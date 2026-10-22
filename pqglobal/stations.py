"""Station selection from FDSN metadata.

Rule (user decision, 2026-10-18; the same rule as the SeisComP San Ramon study): a station is used
when the network itself picked one of its sensors for the catalog's events, that sensor lies within
`radius_km` of the catalog centroid, and it has three components at `min_sps` or more (any sensor
type: broadband, short-period, borehole, accelerometer). One sensor per station: the picked sensor
first in PRECEDENCE, ties broken by the number of picks. Vertical-only sensors are left out (the
picker needs three components), and so is every unpicked sensor.

Writes, into the cluster's station_table/:
  FDSN_station.csv   Network, Code, Latitude, Longitude, Elevation (m)  (the pipeline's table)
  sensors.csv        every picked sensor in the radius with its verdict (used / why not)
  inventory.xml      response-level StationXML of the used sensors (the waveform fetcher rotates
                     with its orientations and divides by its signed sensitivities)
"""
from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
from obspy import Inventory, UTCDateTime
from obspy.geodetics import gps2dist_azimuth

# One sensor per station, first match wins: broadband, short-period, borehole short-period /
# high-rate, low-gain seismometer, accelerometers. Bands not listed come after, velocity first.
PRECEDENCE = ("HH", "EH", "EP", "DP", "HL", "HN", "CN", "DN", "EN", "HG")
BULK_CHUNK = 40           # NCEDC answered 40-sensor bulk requests in ~47 s; 464 at once timed out
BULK_WORKERS = 2          # at most two requests at a time to one data centre


def _rank(band: str) -> tuple:
    if band in PRECEDENCE:
        return (0, PRECEDENCE.index(band))
    return (1, 1 if band[1:2] in ("N", "G") else 0, band)


def is_accelerometer(band: str) -> bool:
    return band[1:2] in ("N", "G")


def _bulk(client, lines, level):
    chunks = [lines[i:i + BULK_CHUNK] for i in range(0, len(lines), BULK_CHUNK)]

    def one(ch):
        for attempt in range(3):
            try:
                return client.get_stations_bulk(ch, level=level)
            except Exception as e:  # noqa: BLE001
                if "No data" in str(e) or "204" in str(e):
                    return Inventory(networks=[], source="")
                if attempt == 2:
                    raise
    inv = Inventory(networks=[], source="")
    with ThreadPoolExecutor(BULK_WORKERS) as ex:
        for part in ex.map(one, chunks):
            inv += part
    return inv


def select_picked(picks_csv: str, centroid, starttime, endtime, out_dir: str,
                  service: str = "NCEDC", radius_km: float = 60.0, min_sps: float = 100.0) -> dict:
    """Apply the station rule; write FDSN_station.csv, sensors.csv, inventory.xml to out_dir."""
    from obspy.clients.fdsn import Client
    client = Client(service, timeout=600)
    t0, t1 = UTCDateTime(starttime), UTCDateTime(endtime)
    lat0, lon0 = centroid

    p = pd.read_csv(picks_csv, dtype={"location": str}, keep_default_na=False)
    p["band"] = p["channel"].str[:2]
    picked = (p.groupby(["network", "station", "location", "band"]).size()
              .rename("n_picks").reset_index())

    # 1. coordinates of every picked station (cheap), then keep those within the radius
    st_inv = _bulk(client, [(n, s, "*", "*", t0, t1)        # (UTCDateTime is unhashable)
                            for n, s in sorted(set(zip(picked.network, picked.station)))],
                   "station")
    dist = {}
    for net in st_inv:
        for sta in net:
            dist[(net.code, sta.code)] = gps2dist_azimuth(lat0, lon0, sta.latitude,
                                                          sta.longitude)[0] / 1e3
    picked["dist_km"] = [dist.get((n, s)) for n, s in zip(picked.network, picked.station)]
    near = picked[picked.dist_km.notna() & (picked.dist_km <= radius_km)].copy()

    # 2. channel + response metadata of the picked sensors within the radius
    lines = [(r.network, r.station, r.location or "--", f"{r.band}?", t0, t1)
             for r in near.itertuples()]
    inv = _bulk(client, lines, "response")

    rows = []
    for r in near.itertuples():
        chans = inv.select(network=r.network, station=r.station, location=r.location,
                           channel=f"{r.band}?", starttime=t0, endtime=t1)
        chs = [c for n in chans for s in n for c in s]
        fast = [c for c in chs if c.sample_rate >= min_sps]
        comps = sorted({c.code[2] for c in fast})
        oriented = all(c.azimuth is not None and c.dip is not None for c in fast)
        z = next((c for c in fast if c.code[2] in ("Z", "1") and c.dip is not None
                  and abs(c.dip) > 60), fast[0] if fast else None)
        if not chs:
            verdict = "no metadata"
        elif not fast:
            verdict = f"below {min_sps:g} sps"
        elif len(comps) < 3:
            verdict = "vertical-only" if len(comps) == 1 else "two components"
        elif not oriented:
            verdict = "no orientation"
        else:
            verdict = "candidate"
        rows.append(dict(network=r.network, station=r.station, location=r.location, band=r.band,
                         sample_rate=max((c.sample_rate for c in chs), default=None),
                         components="".join(comps), dist_km=round(r.dist_km, 3),
                         n_picks=int(r.n_picks), latitude=(z.latitude if z else None),
                         longitude=(z.longitude if z else None),
                         elevation_m=(z.elevation if z else None),
                         kind=("acceleration" if is_accelerometer(r.band) else "velocity"),
                         verdict=verdict))
    sens = pd.DataFrame(rows)
    cand = sens[sens.verdict == "candidate"].copy()
    cand["_rank"] = cand.band.map(_rank)
    cand = cand.sort_values(["network", "station", "_rank", "n_picks"],
                            ascending=[True, True, True, False], kind="mergesort")
    first = ~cand.duplicated(["network", "station"], keep="first")
    sens.loc[cand.index[first], "verdict"] = "used"
    sens.loc[cand.index[~first], "verdict"] = "other sensor of a used station"
    used = sens[sens.verdict == "used"].sort_values(["dist_km"]).reset_index(drop=True)

    os.makedirs(out_dir, exist_ok=True)
    sens.sort_values(["dist_km", "network", "station", "band"]).to_csv(
        os.path.join(out_dir, "sensors.csv"), index=False)
    pd.DataFrame({"Network": used.network, "Code": used.station, "Latitude": used.latitude,
                  "Longitude": used.longitude, "Elevation": used.elevation_m}).to_csv(
        os.path.join(out_dir, "FDSN_station.csv"), index=False)
    keep = Inventory(networks=[], source=inv.source)
    for r in used.itertuples():
        keep += inv.select(network=r.network, station=r.station, location=r.location,
                           channel=f"{r.band}?")
    keep.write(os.path.join(out_dir, "inventory.xml"), format="STATIONXML")
    return dict(picked_stations=int(picked[["network", "station"]].drop_duplicates().shape[0]),
                picked_sensors_in_radius=len(sens),
                used_stations=len(used),
                verdicts=sens.verdict.value_counts().to_dict(),
                used_by_band=used.band.value_counts().to_dict())
