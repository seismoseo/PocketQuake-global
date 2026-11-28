"""The README header figure from the Pinnacles 2026 benchmark run (PyGMT; run after
`examples/pinnacles2026.sh apply`):

  left    the 42 stations, by the model HYPOINVERSE uses for them: NCSN's Pacific-side Gabilan model
          (GAB), or Diablo (DIA) for the stations of NCSN's alternate list and the newer east-side
          stations; the San Andreas fault (USGS NSHM23 traces) and the sequence
  middle  NCSN's catalog (ComCat) against PocketQuake-global's relocation (HypoDD with cross-correlation)
  right   both along strike against depth, from the mainshock

    python examples/pinnacles2026/figure.py [--out docs/images/pinnacles2026.png]
"""
from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

_tmp = tempfile.mkdtemp(prefix="gmt_")                 # ~/.gmt may be read-only
os.environ.setdefault("GMT_SESSIONDIR", _tmp)
import numpy as np
import pandas as pd
import pygmt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
KCR = REPO / "external" / "PocketQuake" / "external" / "korea-cluster-relocation"
STRIKE = 137.12                                        # NC moment tensor, the San Andreas fault here
MAINSHOCK = "nc75468702"
RELOC_COLS = ["id", "lat", "lon", "depth", "x", "y", "z", "ex", "ey", "ez", "yr", "mo", "dy", "hr", "mi",
              "sc", "mag", "nccp", "nccs", "nctp", "ncts", "rcc", "rct", "cid"]


def along(lat, lon, lat0, lon0):
    x = (np.asarray(lon) - lon0) * 111.195 * np.cos(np.radians(lat0))
    y = (np.asarray(lat) - lat0) * 111.195
    s = np.radians(STRIKE)
    return x * np.sin(s) + y * np.cos(s)


def relocated(run_dir, cat):
    """HypoDD dt.cc positions, matched to the ComCat events by origin time (within 3 s)."""
    r = pd.read_csv(run_dir / "2.HypoDD" / "02.dt.cc" / "hypoDD.reloc", sep=r"\s+", header=None, names=RELOC_COLS)
    t = pd.to_datetime(dict(year=r.yr, month=r.mo, day=r.dy, hour=r.hr, minute=r.mi)) + pd.to_timedelta(r.sc, "s")
    ct = pd.to_datetime(cat.time).dt.tz_localize(None)
    k = np.abs(ct.values[None, :] - t.values[:, None]).argmin(1)
    ok = np.abs(ct.values[k] - t.values) < np.timedelta64(3, "s")
    out = r.loc[ok, ["lat", "lon", "depth"]].reset_index(drop=True)
    out["id"], out["mag"] = cat.id.values[k[ok]], cat.mag.values[k[ok]]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", default=str(KCR / "pipeline" / "runs" / "pinnacles2026"))
    ap.add_argument("--cluster-dir", default=str(KCR / "Pinnacles2026_cluster"))
    ap.add_argument("--work", default=str(REPO / "work" / "pinnacles2026"))
    ap.add_argument("--out", default=str(REPO / "docs" / "images" / "pinnacles2026.png"))
    a = ap.parse_args()
    cat = pd.read_csv(HERE / "comcat_sequence.csv")
    rel = relocated(Path(a.run_dir), cat)
    ms = cat[cat.id == MAINSHOCK].iloc[0]
    ms_r = rel[rel.id == MAINSHOCK].iloc[0]
    sta = pd.read_csv(Path(a.cluster_dir) / "station_table" / "FDSN_station.csv")
    alt = {tuple(l.split()[:2]) for l in open(Path(a.work) / "ncsn_alternate_stations.del") if l.strip()}
    sta["dia"] = [(c, n) in alt for c, n in zip(sta.Code, sta.Network)]
    trace = json.load(open(HERE / "usgs_nshm23_san_andreas.geojson"))

    def faults(fig):
        for f in trace["features"]:
            g = f["geometry"]
            for ln in (g["coordinates"] if g["type"] == "MultiLineString" else [g["coordinates"]]):
                ln = np.array(ln)
                fig.plot(x=ln[:, 0], y=ln[:, 1], pen="1.1p,brown")

    size = lambda m: 0.03 + 0.035 * np.clip(m, 0, 5)
    H = 4.6                                            # panel height (cm); map widths follow their regions

    def mercator_width(reg):
        y = lambda lat: np.log(np.tan(np.pi / 4 + np.radians(lat) / 2))
        return H * np.radians(reg[1] - reg[0]) / (y(reg[3]) - y(reg[2]))

    fig = pygmt.Figure()
    title = lambda t: fig.text(text=t, position="TC", justify="BC", offset="0c/0.18c", no_clip=True,
                               font="9p,Helvetica")
    pygmt.config(FONT="8p,Helvetica", FONT_TITLE="9p,Helvetica", MAP_FRAME_TYPE="plain",
                 FORMAT_GEO_MAP="ddd.x", MAP_TITLE_OFFSET="4p")
    blue, orange, red, grey = "0/114/178", "230/159/0", "213/94/0", "gray55"
    # stations
    reg = [-121.95, -120.45, 36.1, 37.15]
    w1 = mercator_width(reg)
    fig.basemap(region=reg, projection=f"M{w1:.3f}c", frame=["WSne", "xa0.5", "ya0.5"])
    fig.coast(shorelines="0.4p,gray50", water="230/238/245", land="white")
    fig.basemap(frame=["WSne", "xa0.5", "ya0.5"], map_scale="jBL+w20k+o0.25c/0.3c+c36.6")   # over the coast
    title("Stations and velocity models")
    faults(fig)
    for dia, col in ((False, blue), (True, orange)):
        s = sta[sta.dia == dia]
        fig.plot(x=s.Longitude, y=s.Latitude, style="i0.24c", fill=col, pen="0.3p,black")
    fig.plot(x=[ms.longitude], y=[ms.latitude], style="a0.42c", fill=red, pen="0.4p,black")
    # NCSN vs PocketQuake-global, map
    fig.shift_origin(xshift=f"{w1 + 1.3:.3f}c")
    reg = [-121.30, -121.12, 36.535, 36.645]
    w2 = mercator_width(reg)
    with pygmt.config(FORMAT_GEO_MAP="ddd.xx"):
        fig.basemap(region=reg, projection=f"M{w2:.3f}c", frame=["WSne", "xa0.1f0.05", "ya0.05f0.025"],
                    map_scale="jBL+w2k+o0.3c/0.45c+c36.6")
    title("NCSN catalog and PocketQuake-global")
    faults(fig)
    fig.plot(x=cat.longitude, y=cat.latitude, style="c", size=size(cat.mag), pen="0.5p," + grey)
    fig.plot(x=rel.lon, y=rel.lat, style="c", size=size(rel.mag), fill=blue, pen="0.2p,black")
    fig.plot(x=[ms_r.lon], y=[ms_r.lat], style="a0.42c", fill=red, pen="0.4p,black")
    # along strike against depth
    fig.shift_origin(xshift=f"{w2 + 1.55:.3f}c")
    w3 = 5.6
    a_cat = along(cat.latitude, cat.longitude, ms.latitude, ms.longitude)
    a_rel = along(rel.lat, rel.lon, ms_r.lat, ms_r.lon)
    fig.basemap(region=[-3, 6, 2, 10], projection=f"X{w3}c/-{H}c",
                frame=["WSne", "xa2f1+lAlong strike from the mainshock (km)", "ya2f1+lDepth (km)"])
    title("Along the fault")
    fig.plot(x=a_cat, y=cat.depth, style="c", size=size(cat.mag), pen="0.5p," + grey)
    fig.plot(x=a_rel, y=rel.depth, style="c", size=size(rel.mag), fill=blue, pen="0.2p,black")
    fig.plot(x=[0], y=[ms_r.depth], style="a0.42c", fill=red, pen="0.4p,black")
    fig.text(x=-2.8, y=9.7, text="NW", font="7p,Helvetica,gray30", justify="LB")
    fig.text(x=5.8, y=9.7, text="SE", font="7p,Helvetica,gray30", justify="RB")
    # one symbol key under the three panels
    total = w1 + 1.3 + w2 + 1.55 + w3
    key = Path(_tmp) / "key.txt"
    key.write_text("\n".join([
        "N 3.5 3.2 2.7 2.4 2.9 2.9",                     # relative column widths, by label length
        f"S 0.2c i 0.24c {blue} 0.3p,black 0.4c Gabilan model (Pacific side)",
        f"S 0.2c i 0.24c {orange} 0.3p,black 0.4c Diablo model (east side)",
        "S 0.2c - 0.4c - 1.1p,brown 0.45c San Andreas fault",
        f"S 0.2c c 0.16c - 0.5p,{grey} 0.4c NCSN catalog",
        f"S 0.2c c 0.16c {blue} 0.2p,black 0.4c PocketQuake-global",
        f"S 0.2c a 0.32c {red} 0.4p,black 0.4c Mw 4.98 mainshock"]) + "\n")
    fig.shift_origin(xshift=f"-{w1 + 1.3 + w2 + 1.55:.3f}c")
    fig.legend(spec=str(key), position=f"JBL+jTL+o0c/1.05c+w{total:.2f}c", box=False)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, dpi=250)
    print(a.out, f"{len(rel)} relocated of {len(cat)}; {int(sta.dia.sum())} of {len(sta)} stations on Diablo")


if __name__ == "__main__":
    main()
