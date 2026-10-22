"""The seismic network's own (analyst) picks for the catalog's events, from an FDSN event service.

Only the picks the network used for its preferred location (its arrivals) count; picks of other
solutions of the same event (e.g. an automatic origin) are left out.

Used for two things:
  * the station rule: the stations the network itself picked for these events (stations.py);
  * an independent check of our automatic picks.

One request per UTC day (the catalog's box, includearrivals=true), cached per day, so a rerun
fetches only the missing days. Each pick is matched to a catalog event by its id (the number
after the network prefix, e.g. nc72540765 = NCEDC event 72540765); events whose ids do not match
fall back to the nearest origin within MATCH_DT_S seconds and MATCH_KM kilometres.

Output columns: event_id (catalog source_id), network, station, location, channel, phase, time
(ISO UTC), uncertainty_s, onset, polarity, evaluation_mode, time_weight, residual_s, distance_km.
"""
from __future__ import annotations

import os
import re
import time as _time

import pandas as pd
from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth

# FDSN event services that return picks (includearrivals=true), by catalog id prefix. USGS ComCat
# itself serves no picks over fdsnws.
PICK_SERVICES = {"nc": "NCEDC", "bk": "NCEDC", "ci": "SCEDC"}
MATCH_DT_S = 2.0
MATCH_KM = 10.0
PAD_DEG = 0.05            # box padding around the catalog (routine vs catalog locations differ)
RETRIES = 3

COLUMNS = ["event_id", "network", "station", "location", "channel", "phase", "time",
           "uncertainty_s", "onset", "polarity", "evaluation_mode", "time_weight", "residual_s",
           "distance_km"]


def service_for(catalog: pd.DataFrame) -> str | None:
    """The pick service for a catalog: the most common id prefix's service, or None."""
    pref = catalog["source_id"].astype(str).str.extract(r"^([a-z]+)", expand=False).dropna()
    return PICK_SERVICES.get(pref.mode().iloc[0]) if len(pref) else None


def _num(eid: str) -> str:
    """'nc72540765' / 'quakeml:nc.anss.org/Event/NC/72540765' -> '72540765'."""
    m = re.search(r"(\d+)$", str(eid))
    return m.group(1) if m else str(eid)


def _day_picks(client, day: UTCDateTime, box) -> pd.DataFrame:
    from obspy.clients.fdsn.header import FDSNNoDataException
    lat0, lat1, lon0, lon1 = box
    for attempt in range(RETRIES):
        try:
            cat = client.get_events(starttime=day, endtime=day + 86400, minlatitude=lat0,
                                    maxlatitude=lat1, minlongitude=lon0, maxlongitude=lon1,
                                    includearrivals=True)
            break
        except FDSNNoDataException:
            return pd.DataFrame(columns=COLUMNS + ["_num", "_otime", "_olat", "_olon"])
        except Exception as e:  # noqa: BLE001  (service hiccups: retry, then give up loudly)
            if attempt == RETRIES - 1:
                raise RuntimeError(f"picks for {day.date}: {e}") from e
            _time.sleep(10 * (attempt + 1))
    rows = []
    for e in cat:
        o = e.preferred_origin() or (e.origins[0] if e.origins else None)
        if o is None:
            continue
        arr = {a.pick_id.id: a for a in o.arrivals}
        for p in e.picks:
            a = arr.get(p.resource_id.id)
            if a is None:              # a pick of another (e.g. automatic) solution of the event
                continue
            w = p.waveform_id
            rows.append(dict(
                event_id="", network=w.network_code, station=w.station_code,
                location=w.location_code or "", channel=w.channel_code,
                phase=(a.phase if a is not None and a.phase else p.phase_hint),
                time=str(p.time), uncertainty_s=p.time_errors.uncertainty, onset=p.onset,
                polarity=p.polarity, evaluation_mode=p.evaluation_mode,
                time_weight=(a.time_weight if a is not None else None),
                residual_s=(a.time_residual if a is not None else None),
                distance_km=(a.distance * 111.19 if a is not None and a.distance is not None
                             else None),
                _num=_num(e.resource_id.id), _otime=str(o.time), _olat=o.latitude,
                _olon=o.longitude))
    return pd.DataFrame(rows, columns=COLUMNS + ["_num", "_otime", "_olat", "_olon"])


def fetch(catalog: pd.DataFrame, out_csv: str, cache_dir: str, service: str | None = None,
          pad_deg: float = PAD_DEG) -> dict:
    """Fetch the network picks of every catalog event; write out_csv. Returns counts."""
    from obspy.clients.fdsn import Client
    service = service or service_for(catalog)
    if service is None:
        raise ValueError("no pick service known for this catalog's ids; pass service=")
    client = Client(service, timeout=600)
    os.makedirs(cache_dir, exist_ok=True)
    times = [UTCDateTime(t) for t in catalog["time"]]
    box = (catalog.latitude.min() - pad_deg, catalog.latitude.max() + pad_deg,
           catalog.longitude.min() - pad_deg, catalog.longitude.max() + pad_deg)
    days = [UTCDateTime(d) for d in sorted({str(t.date) for t in times})]   # UTCDateTime is unhashable
    parts = []
    for i, day in enumerate(days, 1):
        f = os.path.join(cache_dir, f"picks_{day.strftime('%Y%m%d')}.csv")
        if not os.path.exists(f):
            t0 = _time.time()
            _day_picks(client, day, box).to_csv(f + ".part", index=False)
            os.replace(f + ".part", f)
            print(f"[picks] {day.date} ({i}/{len(days)}) {_time.time() - t0:.0f} s", flush=True)
        parts.append(pd.read_csv(f, dtype={"location": str, "_num": str}, keep_default_na=False,
                                 na_values={"uncertainty_s": [""], "time_weight": [""],
                                            "residual_s": [""], "distance_km": [""]}))
    allp = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=COLUMNS)

    # match service events to catalog events: by id number, else nearest origin in time and space
    cat = catalog.assign(_num=catalog["source_id"].map(_num),
                         _t=[UTCDateTime(t) for t in catalog["time"]])
    by_num = dict(zip(cat["_num"], cat["source_id"]))
    svc = allp.drop_duplicates("_num")[["_num", "_otime", "_olat", "_olon"]]
    match, n_time = {}, 0
    for num, otime, olat, olon in zip(svc["_num"], svc["_otime"], svc["_olat"], svc["_olon"]):
        if num in by_num:
            match[num] = by_num[num]
            continue
        t = UTCDateTime(otime)
        near = cat[[abs(x - t) <= MATCH_DT_S for x in cat["_t"]]]
        if len(near):
            d = [gps2dist_azimuth(olat, olon, la, lo)[0] / 1e3
                 for la, lo in zip(near.latitude, near.longitude)]
            j = min(range(len(d)), key=d.__getitem__)
            if d[j] <= MATCH_KM:
                match[num] = near.iloc[j]["source_id"]
                n_time += 1
    allp["event_id"] = allp["_num"].map(match)
    picks = allp[allp["event_id"].notna()][COLUMNS]
    os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
    picks.to_csv(out_csv, index=False)
    n_ev = picks["event_id"].nunique()
    return dict(service=service, days=len(days), picks=len(picks), events_with_picks=n_ev,
                catalog_events=len(catalog), matched_by_time=n_time,
                service_events_unmatched=int(svc.shape[0] - len(match)))
