"""The benchmark numbers of a Pinnacles 2026 run (after `examples/pinnacles2026.sh apply`), beside the
reference values of PocketQuake-global v0.3.0:

  picks       PhaseNet+ against NCSN's analyst picks at the same station and phase (event service picks with
              evaluation mode 'manual'): median |difference| and the share within 0.1 s
  locations   events located by HYPOINVERSE and relocated by HypoDD (dt.ct, dt.cc); HYPOINVERSE median RMS
  mechanism   SKHASH: the mainshock's Kagan angle to the NC moment tensor (137.12/72.56/175.92); events whose
              (first) solution has quality A or B, and those of them with a single solution
  shape       HypoDD dt.cc, and NCSN's catalog for comparison: the line through the epicenters within 6 km of
              the mainshock (principal axis): its strike, the scatter across it (sd), the mainshock's distance
              from it; aftershocks south-east of the mainshock along strike 137.12 deg; distance south-west
              of the USGS NSHM23 fault trace
  ldeo_ddrt   against LDEO's real-time double-difference catalog (an independent relocation; the excerpt
              ldeo_ddrt_sequence.txt), on the events both relocated: the offset between the two clouds, the
              differences left after removing it, the width of the line in each (compare_ddrt.py), the plane
              fitted to the aftershocks and its angle to the USGS Mww and NC moment-tensor planes

    python examples/pinnacles2026/check.py
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # the repository, when run as a script
import pqglobal                                        # noqa: F401,E402  (sets the PocketQuake paths)
from pqglobal.waveforms import event_id

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
KCR = REPO / "external" / "PocketQuake" / "external" / "korea-cluster-relocation"
MAINSHOCK = "nc75468702"
NC_MT = (137.12, 72.56, 175.92)
STRIKE = 137.12
RELOC_COLS = ["id", "lat", "lon", "depth", "x", "y", "z", "ex", "ey", "ez", "yr", "mo", "dy", "hr", "mi",
              "sc", "mag", "nccp", "nccs", "nctp", "ncts", "rcc", "rct", "cid"]
# v0.3.0, one run from an empty folder on 2026-11-28 (PyTorch 2.11 with CUDA 12.8, NVIDIA RTX PRO 6000 Blackwell)
REFERENCE = json.load(open(HERE / "reference.json")) if (HERE / "reference.json").is_file() else {}


def km(lat, lon, lat0, lon0):
    return (np.asarray(lon) - lon0) * 111.195 * np.cos(np.radians(lat0)), (np.asarray(lat) - lat0) * 111.195


def kagan(sdr1, sdr2):
    """Smallest rotation (deg) between two double couples (strike, dip, rake)."""
    def frame(s, d, r):
        s, d, r = np.radians([s, d, r])
        n = np.array([-np.sin(d) * np.sin(s), np.sin(d) * np.cos(s), -np.cos(d)])
        u = np.array([np.cos(r) * np.cos(s) + np.cos(d) * np.sin(r) * np.sin(s),
                      np.cos(r) * np.sin(s) - np.cos(d) * np.sin(r) * np.cos(s), -np.sin(d) * np.sin(r)])
        t, p = (n + u) / np.sqrt(2), (n - u) / np.sqrt(2)
        return np.c_[t, p, np.cross(t, p)]
    A, B = frame(*sdr1), frame(*sdr2)
    return min(float(np.degrees(np.arccos(np.clip((np.trace(A @ Q @ B.T) - 1) / 2, -1, 1))))
               for Q in (np.diag([1, 1, 1]), np.diag([1, -1, -1]), np.diag([-1, 1, -1]), np.diag([-1, -1, 1])))


def trace_offset(lat, lon, lat0, lon0):
    """Signed distance (km, + = south-west) from the nearest segment of the NSHM23 San Andreas trace."""
    g = json.load(open(HERE / "usgs_nshm23_san_andreas.geojson"))
    segs = []
    for f in g["features"]:
        geo = f["geometry"]
        for ln in (geo["coordinates"] if geo["type"] == "MultiLineString" else [geo["coordinates"]]):
            ln = np.array(ln)
            x, y = km(ln[:, 1], ln[:, 0], lat0, lon0)
            segs += [(x[i], y[i], x[i + 1], y[i + 1]) for i in range(len(x) - 1)]
    out = []
    for X, Y in zip(*km(lat, lon, lat0, lon0)):
        best = None
        for x1, y1, x2, y2 in segs:
            dx, dy = x2 - x1, y2 - y1
            t = np.clip(((X - x1) * dx + (Y - y1) * dy) / (dx * dx + dy * dy), 0, 1)
            d = np.hypot(X - (x1 + t * dx), Y - (y1 + t * dy))
            if best is None or d < best[0]:
                best = (d, np.sign(dx * (Y - y1) - dy * (X - x1)), dy)
        d, side, dy = best
        out.append(d * (-side if dy < 0 else side))
    return np.array(out)


def reloc(path, cat):
    r = pd.read_csv(path, sep=r"\s+", header=None, names=RELOC_COLS)
    t = pd.to_datetime(dict(year=r.yr, month=r.mo, day=r.dy, hour=r.hr, minute=r.mi)) + pd.to_timedelta(r.sc, "s")
    ct = pd.to_datetime(cat.time).dt.tz_localize(None)
    k = np.abs(ct.values[None, :] - t.values[:, None]).argmin(1)
    ok = np.abs(ct.values[k] - t.values) < np.timedelta64(3, "s")
    r = r[ok].reset_index(drop=True)
    r["event"], r["time"] = cat.source_id.values[k[ok]], ct.values[k[ok]]
    return r


def shape(lat, lon, ids, times):
    """The line through the epicenters within 6 km of the mainshock (principal axis): its strike, the scatter
    across it, the mainshock's distance from it (+ = south-west); aftershocks south-east of the mainshock
    along strike 137.12 deg; median distance south-west of the fault trace."""
    h = np.flatnonzero(ids == MAINSHOCK)[0]
    x, y = km(lat, lon, lat[h], lon[h])
    near = (np.hypot(x, y) < 6.0) & (ids != MAINSHOCK)
    P = np.c_[x[near], y[near]]
    c = P.mean(0)
    v = np.linalg.svd(P - c)[2][0]
    nrm = np.array([v[1], -v[0]])
    nrm = nrm if nrm @ np.array([-np.sin(np.radians(47.12)), -np.cos(np.radians(47.12))]) > 0 else -nrm
    st = np.radians(STRIKE)
    along = x * np.sin(st) + y * np.cos(st)
    after = times > times[h]
    return dict(line_strike_deg=round(float(np.degrees(np.arctan2(v[0], v[1])) % 180), 1),
                across_line_sd_m=round(float(((P - c) @ nrm).std() * 1e3)),
                mainshock_off_line_m=round(float((-c) @ nrm * 1e3)),
                aftershocks_southeast=round(float((along[after] > 0).mean()), 3),
                km_southwest_of_trace_median=round(float(np.median(trace_offset(lat, lon, lat[h], lon[h]))), 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", default=str(KCR / "pipeline" / "runs" / "pinnacles2026"))
    ap.add_argument("--work", default=str(REPO / "work" / "pinnacles2026"))
    ap.add_argument("--write-reference", action="store_true", help="store these values as the reference")
    a = ap.parse_args()
    run, work = Path(a.run_dir), Path(a.work)
    need = [work / "event_catalog.csv", work / "network_picks.csv", run / "picks",
            run / "2.HypoDD" / "01.dt.ct" / "hypoDD.reloc", run / "2.HypoDD" / "02.dt.cc" / "hypoDD.reloc"]
    need += [p for p in (glob.glob(str(run / "1.HypoInv" / "*" / "*.sum"))[:1] or [run / "1.HypoInv/<model>/*.sum"])]
    need += [p for p in (glob.glob(str(run / "3.FocalMech" / "*" / "mechanisms.csv"))[:1] or [run / "3.FocalMech/<model>/mechanisms.csv"])]
    missing = [str(p) for p in need if not Path(p).exists()]
    if missing:
        raise SystemExit("the run is not finished (run examples/pinnacles2026.sh apply first); missing:\n  "
                         + "\n  ".join(missing))
    if not REFERENCE and not a.write_reference:
        print(f"note: no {HERE / 'reference.json'}; the reference column stays empty")
    cat = pd.read_csv(work / "event_catalog.csv")
    cat["eid"] = [event_id(t) for t in cat.time]
    res = {}
    # picks
    p = pd.concat([pd.read_csv(f, dtype={"Event_ID": str, "Station": str}) for f in glob.glob(str(run / "picks" / "*_picks.csv"))])
    p = p.merge(cat[["eid", "source_id"]], left_on="Event_ID", right_on="eid")
    n = pd.read_csv(work / "network_picks.csv", dtype={"station": str})
    n = n[n.evaluation_mode == "manual"]
    n = n.assign(t=pd.to_datetime(n.time)).groupby(["event_id", "network", "station", "phase"]).t.min().reset_index()
    m = p.assign(t=pd.to_datetime(p.Time)).merge(n, left_on=["source_id", "Network", "Station", "Phase"],
                                                  right_on=["event_id", "network", "station", "phase"], suffixes=("", "_ncsn"))
    d = (m.t - m.t_ncsn).dt.total_seconds().abs().round(3)     # picks are on a 10-ms grid: no float edge at 0.1 s
    for ph in ("P", "S"):
        x = d[m.Phase == ph]
        res[f"picks_{ph}"] = dict(pairs=int(len(x)), median_abs_s=round(float(x.median()), 3),
                                  within_0p1=round(float((x <= 0.1).mean()), 3))
    # locations
    from pipeline.core import sumio
    s = sumio.read_sum(glob.glob(str(run / "1.HypoInv" / "*" / "*.sum"))[0])
    ct = reloc(run / "2.HypoDD" / "01.dt.ct" / "hypoDD.reloc", cat)
    cc = reloc(run / "2.HypoDD" / "02.dt.cc" / "hypoDD.reloc", cat)
    res["located"] = dict(events=int(len(cat)), hypoinverse=int(len(s)), hypodd_dtct=int(len(ct)), hypodd_dtcc=int(len(cc)),
                          hypoinverse_median_rms_s=round(float(s.rms.median()), 3))
    # mechanisms
    fm = pd.read_csv(glob.glob(str(run / "3.FocalMech" / "*" / "mechanisms.csv"))[0], dtype={"event_id": str})
    ms_fm = fm[fm.event_id == cat.set_index("source_id").loc[MAINSHOCK, "eid"]].iloc[0]
    first = fm.drop_duplicates("event_id")              # SKHASH lists further solutions of an event below its first
    single = first[~first.mult_solution_flag.astype(bool)]
    res["mechanisms"] = dict(events=int(len(first)), quality_A_or_B=int(first.quality.isin(["A", "B"]).sum()),
                             quality_A_or_B_single_solution=int(single.quality.isin(["A", "B"]).sum()),
                             mainshock=[float(ms_fm.strike), float(ms_fm.dip), float(ms_fm.rake)],
                             mainshock_kagan_to_nc_mt_deg=round(kagan((ms_fm.strike, ms_fm.dip, ms_fm.rake), NC_MT), 1))
    # shape of the dt.cc relocation, and of NCSN's catalog (ComCat) for comparison
    res["shape_dtcc"] = shape(cc.lat.values, cc.lon.values, cc.event.values, cc.time.values)
    res["shape_ncsn_catalog"] = shape(cat.latitude.values, cat.longitude.values, cat.source_id.values,
                                      pd.to_datetime(cat.time).dt.tz_localize(None).values)
    # against LDEO's real-time double-difference catalog
    import compare_ddrt as cd
    r = cd.compare(cc[["event", "lat", "lon", "depth"]], pd.read_csv(HERE / "comcat_sequence.csv"),
                   cd.read_ddrt(HERE / "ldeo_ddrt_sequence.txt"), cd.trace_segments(HERE / "usgs_nshm23_san_andreas.geojson"))[0]
    o, d, rel = r["catalogs"]["ours"], r["catalogs"]["ddrt"], r["relative_after_shift"]
    res["ldeo_ddrt"] = dict(compared=r["compared"], offset_horizontal_m=r["shift_ddrt_minus_ours_m"]["horizontal"],
                            offset_down_m=r["shift_ddrt_minus_ours_m"]["down"],
                            along_median_m=rel["along_strike"]["median_m"], across_median_m=rel["across_strike"]["median_m"],
                            depth_median_m=rel["depth"]["median_m"], width_m=[o["width_m"], d["width_m"]],
                            off_line_gt_0p5km=[o["off_line_gt_0p5km"], d["off_line_gt_0p5km"]],
                            ldeo_off_line_from_preliminary_input=d["off_line_from_preliminary_input"],
                            aftershocks=r["aftershocks"], plane=[o["plane_strike"], o["plane_dip"]],
                            plane_dip_95pct=o["plane_dip_95pct"],
                            plane_to_usgs_mww_deg=o["angle_to_moment_tensors"]["USGS Mww"]["deg"],
                            plane_to_nc_mt_deg=o["angle_to_moment_tensors"]["NC moment tensor"]["deg"])
    if a.write_reference:
        json.dump(res, open(HERE / "reference.json", "w"), indent=1)
    rows = []
    for k, sub in res.items():
        for q, val in sub.items():
            ref = REFERENCE.get(k, {}).get(q, "")
            rows.append((f"{k}.{q}", val, ref))
    w = max(len(r[0]) for r in rows)
    print(f"{'quantity':{w}}  {'this run':>22}  {'v0.3.0 reference':>22}")
    for q, val, ref in rows:
        print(f"{q:{w}}  {str(val):>22}  {str(ref):>22}")


if __name__ == "__main__":
    main()
