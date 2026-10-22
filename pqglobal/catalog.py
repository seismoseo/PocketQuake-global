"""Event catalogs: a USGS ComCat CSV or an FDSN event query → the cluster's UTC catalog CSV.

Output columns (what the relocation pipeline's UTC reader takes): time (ISO 8601 UTC), latitude,
longitude, depth (km), magnitude, plus mag_type and source_id for the record.

Rule: the pipeline names each event's waveform folder by its UTC second (YYYYMMDDHHMMSS), so two
events in the same UTC second would share one folder. Of such a pair the smaller event is dropped;
every drop is written to <out>.dropped.csv and reported.
"""
from __future__ import annotations

import os

import pandas as pd
from obspy import UTCDateTime

COLUMNS = ["time", "latitude", "longitude", "depth", "magnitude", "mag_type", "source_id"]


def from_usgs_csv(path: str) -> pd.DataFrame:
    """ComCat CSV (earthquake.usgs.gov search 'CSV' export): time, latitude, longitude, depth, mag,
    magType, id, ..."""
    df = pd.read_csv(path)
    need = {"time", "latitude", "longitude", "depth", "mag"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"{path}: not a USGS ComCat CSV (missing {sorted(missing)})")
    return pd.DataFrame({
        "time": [str(UTCDateTime(t)) for t in df["time"]],
        "latitude": df["latitude"].astype(float), "longitude": df["longitude"].astype(float),
        "depth": df["depth"].astype(float), "magnitude": df["mag"].astype(float),
        "mag_type": df.get("magType", pd.Series([""] * len(df))).astype(str),
        "source_id": df.get("id", pd.Series([""] * len(df))).astype(str),
    })


def from_usgs(starttime, endtime, box=None, latitude: float | None = None,
              longitude: float | None = None, maxradius_km: float | None = None,
              minmagnitude: float | None = None, cache_csv: str | None = None) -> pd.DataFrame:
    """USGS ComCat query in its CSV format (QuakeML for hundreds of events takes many minutes)."""
    import urllib.parse
    import urllib.request
    q = dict(starttime=str(UTCDateTime(starttime)), endtime=str(UTCDateTime(endtime)),
             format="csv", orderby="time-asc")
    if box is not None:
        q.update(minlatitude=box[0], maxlatitude=box[1], minlongitude=box[2], maxlongitude=box[3])
    else:
        q.update(latitude=latitude, longitude=longitude, maxradius=maxradius_km / 111.19)
    if minmagnitude is not None:
        q["minmagnitude"] = minmagnitude
    url = "https://earthquake.usgs.gov/fdsnws/event/1/query?" + urllib.parse.urlencode(q)
    with urllib.request.urlopen(url, timeout=600) as r:
        text = r.read().decode()
    if cache_csv:                                    # keep the raw ComCat answer for the record
        with open(cache_csv, "w") as f:
            f.write(text)
    import io
    return from_usgs_csv(io.StringIO(text))


def from_fdsn(client: str, starttime, endtime, latitude: float | None = None,
              longitude: float | None = None, maxradius_km: float | None = None,
              minmagnitude: float | None = None, box=None) -> pd.DataFrame:
    """FDSN event query (e.g. client="USGS") around a point, or in box=(lat0, lat1, lon0, lon1)."""
    from obspy.clients.fdsn import Client
    kw = dict(starttime=UTCDateTime(starttime), endtime=UTCDateTime(endtime), orderby="time-asc")
    if box is not None:
        kw.update(minlatitude=box[0], maxlatitude=box[1], minlongitude=box[2], maxlongitude=box[3])
    else:
        kw.update(latitude=latitude, longitude=longitude, maxradius=maxradius_km / 111.19)
    if minmagnitude is not None:
        kw["minmagnitude"] = minmagnitude
    cat = Client(client, timeout=120).get_events(**kw)
    rows = []
    for e in cat:
        o = e.preferred_origin() or e.origins[0]
        m = e.preferred_magnitude() or (e.magnitudes[0] if e.magnitudes else None)
        rows.append(dict(time=str(o.time), latitude=o.latitude, longitude=o.longitude,
                         depth=(o.depth or 0.0) / 1000.0, magnitude=(m.mag if m else float("nan")),
                         mag_type=(m.magnitude_type if m else ""),
                         source_id=str(e.resource_id).split("eventid=")[-1].split("&")[0]))
    return pd.DataFrame(rows, columns=COLUMNS)


def write(df: pd.DataFrame, out_csv: str) -> dict:
    """Sort by time, apply the same-second rule, write. Returns {'kept': n, 'dropped': n}."""
    df = df.copy()
    df["_t"] = [UTCDateTime(t) for t in df["time"]]
    df["_sec"] = [t.strftime("%Y%m%d%H%M%S") for t in df["_t"]]
    df = df.sort_values(["_sec", "magnitude"], ascending=[True, False], kind="mergesort")
    dup = df.duplicated("_sec", keep="first")
    dropped, kept = df[dup], df[~dup].sort_values("_t", kind="mergesort")
    os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
    kept[COLUMNS].to_csv(out_csv, index=False)
    if len(dropped):
        dropped[COLUMNS].to_csv(out_csv[:-4] + ".dropped.csv", index=False)
    return {"kept": len(kept), "dropped": len(dropped)}


def centroid_and_bounds(df: pd.DataFrame, pad_deg: float = 0.2):
    """Epicenter = catalog centroid; region bounds = catalog bounding box + pad (as pocketquake.sh)."""
    lat, lon = float(df.latitude.mean()), float(df.longitude.mean())
    b = (float(df.latitude.min()) - pad_deg, float(df.latitude.max()) + pad_deg,
         float(df.longitude.min()) - pad_deg, float(df.longitude.max()) + pad_deg)
    return (round(lat, 4), round(lon, 4)), tuple(round(v, 4) for v in b)
