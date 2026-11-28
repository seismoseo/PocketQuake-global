"""PocketQuake-global's relocation of the Pinnacles 2026 sequence against LDEO's real-time double-difference
catalog (DDRT; Waldhauser 2009) and NCSN's routine catalog (ComCat), on the same events.

DDRT relocates each new NCSN event, from NCSN's picks and waveform cross-correlation, relative to the events of
LDEO's 1984-2021 double-difference base catalog (Waldhauser & Schaff 2008). It is computed independently of
PocketQuake-global: other picks, other station set, other velocity model, other reference events.

Numbers (JSON) and one row per event (CSV):
  matched       events in both catalogs; DDRT events without a DD solution (METH = HINV) are left out
  input         DDRT's VER field: the NCSN solution version it started from (0 = NCSN's first, preliminary one)
  shift         median position of DDRT minus ours (east, north, down): where the two clouds sit
  relative      after removing that shift, |DDRT - ours| along strike, across strike and in depth
  width         spread across the line (robust: IQR / 1.349) and events more than 0.5 km off it, events along
                the section stretch of the line, for NCSN, DDRT and ours
  plane         plane fitted to the aftershocks (SVD; strike/dip by the right-hand rule), each catalog with its own
                positions of the same events; dip 95 % range from resampling the aftershocks; angle to the
                mainshock's NW-SE moment-tensor planes (NC Mw: 137.12/72.56; USGS Mww: 135.04/81.98)
Figures (PyGMT; --no-figures skips them): <out>/ldeo_ddrt_comparison.png (map and across-fault section of each catalog)
and <out>/ldeo_ddrt_differences.png (ours and DDRT after removing the shift, tied event by event).

    python examples/pinnacles2026/compare_ddrt.py      # after examples/pinnacles2026.sh apply
Options: --reloc another relocation (a hypoDD.reloc, or a CSV with event, lat, lon, depth and optionally
solution); --classes a CSV with id and class to define the aftershocks (default: after the mainshock and
within 7 km of it); --out-dir (default work/pinnacles2026/ldeo_ddrt). The README shows the reference
run's copy of the comparison figure (docs/images/pinnacles2026_ldeo_ddrt.png).
"""
from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
MAINSHOCK = "nc75468702"
STRIKE = 137.12                                       # NC moment tensor, NW-SE plane: the strike frame
MT_PLANES = {"NC moment tensor": (137.12, 72.56),     # ComCat moment-tensor products of nc75468702
             "USGS Mww": (135.04, 81.98)}
RELOC_COLS = ["id", "lat", "lon", "depth", "x", "y", "z", "ex", "ey", "ez", "yr", "mo", "dy", "hr", "mi",
              "sc", "mag", "nccp", "nccs", "nctp", "ncts", "rcc", "rct", "cid"]
DDRT_COLS = "yr mo dy hr mi sc lat lon dep ex ey ez mag cusp ver base meth".split()
ZONE_KM = 7.0                                         # default aftershock rule: after the mainshock, within 7 km


# ----------------------------------------------------------------------------------------------- inputs
def read_ddrt(path):
    rows = [ln.split() for ln in open(path, errors="replace") if ln[:2] in ("19", "20") and len(ln.split()) == 17]
    d = pd.DataFrame(rows, columns=DDRT_COLS)
    for c in DDRT_COLS[:14]:
        d[c] = pd.to_numeric(d[c])
    d["event"] = "nc" + d.cusp.astype(int).astype(str)
    d["input_version"] = d.ver.str.split(".").str[0].astype(int)
    return d


def read_ours(path, cat, solution="hypodd_cc"):
    """Our positions by event id: a location table with an `event` column (and `solution`), or a hypoDD.reloc
    matched to the catalog by origin time (within 3 s)."""
    path = Path(path)
    if path.suffix == ".csv":
        r = pd.read_csv(path)
        if "solution" in r:
            r = r[r.solution == solution]
        return r[["event", "lat", "lon", "depth"]].reset_index(drop=True)
    r = pd.read_csv(path, sep=r"\s+", header=None, names=RELOC_COLS)
    t = pd.to_datetime(dict(year=r.yr, month=r.mo, day=r.dy, hour=r.hr, minute=r.mi)) + pd.to_timedelta(r.sc, "s")
    ct = pd.to_datetime(cat.time).dt.tz_localize(None)
    k = np.abs(ct.values[None, :] - t.values[:, None]).argmin(1)
    ok = np.abs(ct.values[k] - t.values) < np.timedelta64(3, "s")
    out = r.loc[ok, ["lat", "lon", "depth"]].reset_index(drop=True)
    out["event"] = cat.id.values[k[ok]]
    return out


def trace_segments(path):
    g = json.load(open(path))
    segs = []
    for f in g["features"]:
        geo = f["geometry"]
        for ln in (geo["coordinates"] if geo["type"] == "MultiLineString" else [geo["coordinates"]]):
            segs.append(np.array(ln)[:, :2])
    return segs


# ----------------------------------------------------------------------------------------------- geometry
class Frame:
    """Local km frame at the NCSN mainshock epicenter: east, north, down; along strike (SE +), across (NE +)."""
    def __init__(self, lat0, lon0):
        self.lat0, self.lon0 = lat0, lon0
        s = np.radians(STRIKE)
        self.along_u = np.array([np.sin(s), np.cos(s)])
        self.ne_u = np.array([-np.cos(s), np.sin(s)])

    def en(self, lat, lon):
        return np.c_[(np.asarray(lon) - self.lon0) * 111.195 * np.cos(np.radians(self.lat0)),
                     (np.asarray(lat) - self.lat0) * 111.195]

    def ena(self, lat, lon, dep):
        return np.c_[self.en(lat, lon), np.asarray(dep, float)]


def ne_of_trace(P, segs, fr):
    """Distance (km) north-east of the nearest segment of the fault trace (+ = NE side)."""
    S = []
    for ln in segs:
        xy = fr.en(ln[:, 1], ln[:, 0])
        S += [(xy[i], xy[i + 1]) for i in range(len(xy) - 1)]
    out = []
    for p in P[:, :2]:
        best = None
        for a, b in S:
            v = b - a
            t = np.clip((p - a) @ v / (v @ v), 0, 1)
            d = np.hypot(*(p - a - t * v))
            if best is None or d < best[0]:
                side = np.sign(fr.ne_u @ (p - a - t * v)) or 1.0
                best = (d, side)
        out.append(best[0] * best[1])
    return np.array(out)


def plane_fit(P):
    c = P.mean(0)
    n = np.linalg.svd(P - c)[2][2]
    if n[2] > 0:
        n = -n                                         # upward normal: tilts toward the dip direction
    strike = (np.degrees(np.arctan2(n[0], n[1])) - 90) % 360
    return strike, float(np.degrees(np.arccos(-n[2]))), n, c


def normal(strike, dip):
    """Upward unit normal (east, north, down) of a plane given by the right-hand rule."""
    s, d = np.radians(strike), np.radians(dip)
    return np.array([np.sin(d) * np.cos(s), -np.sin(d) * np.sin(s), -np.cos(d)])


def angle(n1, n2):
    return float(np.degrees(np.arccos(min(1.0, abs(float(n1 @ n2))))))


def robust_width(x):
    x = x - np.median(x)
    return float(np.subtract(*np.percentile(x, [75, 25])) / 1.349), int((np.abs(x) > 0.5).sum())


# ----------------------------------------------------------------------------------------------- comparison
def compare(ours, cat, ddrt, segs, classes=None, section=(-4.0, 7.0), nboot=2000, seed=1128):
    """The numbers (dict) and the per-event table. `ours`: event, lat, lon, depth; `cat`: ComCat CSV
    (id, time, latitude, longitude, depth, mag); `ddrt`: read_ddrt(); `classes`: optional id -> class."""
    cat = cat.rename(columns={"latitude": "ncsn_lat", "longitude": "ncsn_lon", "depth": "ncsn_depth"})
    m = ours.merge(cat[["id", "time", "mag", "ncsn_lat", "ncsn_lon", "ncsn_depth"]], left_on="event", right_on="id")
    m = m.merge(ddrt[["event", "lat", "lon", "dep", "ex", "ey", "ez", "input_version", "meth"]].rename(
        columns={"lat": "dd_lat", "lon": "dd_lon", "dep": "dd_depth"}), on="event", how="left")
    res = dict(ours=int(len(ours)), in_ddrt=int(m.meth.notna().sum()),
               ddrt_without_dd_solution=int((m.meth == "HINV").sum()))
    m = m[m.meth.astype(str).str.startswith("DDRT")].reset_index(drop=True)
    res["compared"] = int(len(m))
    c0 = cat.set_index("id").loc[MAINSHOCK]
    fr = Frame(c0.ncsn_lat, c0.ncsn_lon)
    t_ms = pd.to_datetime(c0.time)
    after = (pd.to_datetime(m.time) > t_ms).values
    tt = pd.to_datetime(m.time)
    res["ddrt_input_version"] = {lab: {str(k): int(v) for k, v in m.loc[sel, "input_version"].value_counts().sort_index().items()}
                                 for lab, sel in (("before_mainshock", tt < t_ms), ("mainshock", m.event == MAINSHOCK),
                                                  ("after_mainshock", tt > t_ms))}
    X = {"ncsn": fr.ena(m.ncsn_lat, m.ncsn_lon, m.ncsn_depth), "ddrt": fr.ena(m.dd_lat, m.dd_lon, m.dd_depth),
         "ours": fr.ena(m.lat, m.lon, m.depth)}
    h = int(np.flatnonzero(m.event.values == MAINSHOCK)[0])
    # where the clouds sit, and the relative positions after removing that
    d = X["ddrt"] - X["ours"]
    sh = np.median(d, 0)
    r = d - sh
    along_d, across_d = r[:, :2] @ fr.along_u, r[:, :2] @ fr.ne_u
    q = lambda v: dict(median_m=round(float(np.median(np.abs(v)) * 1e3)), p90_m=round(float(np.percentile(np.abs(v), 90) * 1e3)))
    res["shift_ddrt_minus_ours_m"] = dict(east=round(float(sh[0] * 1e3)), north=round(float(sh[1] * 1e3)),
                                          down=round(float(sh[2] * 1e3)), horizontal=round(float(np.hypot(*sh[:2]) * 1e3)))
    res["relative_after_shift"] = dict(along_strike=q(along_d), across_strike=q(across_d), depth=q(r[:, 2]))
    res["ddrt_formal_95pct_error_median_m"] = dict(east=round(float(m.ex.median() * 1e3)),
                                                   north=round(float(m.ey.median() * 1e3)),
                                                   down=round(float(m.ez.median() * 1e3)))
    res["mainshock_ddrt_minus_ours_m"] = dict(horizontal=round(float(np.hypot(*d[h, :2]) * 1e3)), down=round(float(d[h, 2] * 1e3)),
                                              horizontal_after_shift=round(float(np.hypot(*r[h, :2]) * 1e3)),
                                              down_after_shift=round(float(r[h, 2] * 1e3)))
    # per catalog: distance from the trace, width of the line, fault plane of the aftershocks
    if classes is not None:
        aft = m.event.map(classes).eq("aftershock").values
    else:
        aft = after & (np.linalg.norm(X["ours"] - X["ours"][h], axis=1) < ZONE_KM)
    rng = np.random.default_rng(seed)
    idx = [rng.integers(0, aft.sum(), aft.sum()) for _ in range(nboot)]
    res["aftershocks"] = int(aft.sum())
    per = {}
    for k, P in X.items():
        a = (P[:, :2] - P[h, :2]) @ fr.along_u
        sec = (a >= section[0]) & (a <= section[1])
        ne = ne_of_trace(P, segs, fr)
        w, off = robust_width(ne[sec])
        st, dp, n, c = plane_fit(P[aft])
        boot = [plane_fit(P[aft][i]) for i in idx]
        dips = np.array([b[1] for b in boot])
        to_mt = {}
        for lab, (s_mt, d_mt) in MT_PLANES.items():        # angle between the planes (normals), and its 95 % range
            nm = normal(s_mt, d_mt)
            to_mt[lab] = dict(deg=round(angle(n, nm), 1),
                              deg_95pct=[round(float(x), 1) for x in np.percentile([angle(b[2], nm) for b in boot], [2.5, 97.5])])
        per[k] = dict(km_ne_of_trace_median=round(float(np.median(ne[sec])), 2), width_m=round(w * 1e3),
                      off_line_gt_0p5km=off, events_in_section=int(sec.sum()),
                      plane_strike=round(float(st), 1), plane_dip=round(dp, 1),
                      plane_dip_95pct=[round(float(x), 1) for x in np.percentile(dips, [2.5, 97.5])],
                      off_plane_sd_m=round(float(np.std((P[aft] - c) @ n) * 1e3)), angle_to_moment_tensors=to_mt)
        if k == "ddrt":                                # off-line events whose DDRT input was NCSN's first solution
            far = np.abs(ne[sec] - np.median(ne[sec])) > 0.5
            per[k]["off_line_from_preliminary_input"] = int((far & (m.input_version.values[sec] == 0)).sum())
        m[f"{k}_ne_of_trace_km"] = ne
        m[f"{k}_along_km"] = (P[:, :2] - X["ours"][h, :2]) @ fr.along_u
    res["catalogs"] = per
    m["aftershock"] = aft
    m["ddrt_minus_ours_after_shift_along_m"] = np.round(along_d * 1e3)
    m["ddrt_minus_ours_after_shift_across_m"] = np.round(across_d * 1e3)
    m["ddrt_minus_ours_after_shift_down_m"] = np.round(r[:, 2] * 1e3)
    return res, m, dict(frame=fr, X=X, h=h, aft=aft, shift=sh, section=section)


# ----------------------------------------------------------------------------------------------- figures
COL = {"ncsn": "153/153/153", "ddrt": "0/158/115", "ours": "0/114/178"}     # Okabe-Ito
NAME = {"ncsn": "NCSN catalog", "ddrt": "LDEO real-time DD", "ours": "PocketQuake-global"}
RED = "213/94/0"


def figures(res, m, g, segs, out_dir):
    tmp = tempfile.mkdtemp(prefix="gmt_")
    os.environ.setdefault("GMT_SESSIONDIR", tmp)
    os.environ.setdefault("GMT_TMPDIR", tmp)
    import pygmt
    fr, X, h, aft = g["frame"], g["X"], g["h"], g["aft"]
    size = lambda mag: 0.03 + 0.03 * np.clip(np.asarray(mag, float), 0, 5)
    mag = m.mag.values
    letter = lambda fig, t, dx: fig.text(text=t, position="TL", justify="BL", offset=f"{dx}c/0.25c", no_clip=True,
                                         font="10p,Helvetica-Bold")
    title = lambda fig, t: fig.text(text=t, position="TC", justify="BC", offset="0c/0.2c", no_clip=True,
                                    font="9p,Helvetica")

    def ll(P):                                          # km frame -> lon, lat
        return (fr.lon0 + P[:, 0] / (111.195 * np.cos(np.radians(fr.lat0))), fr.lat0 + P[:, 1] / 111.195)

    # ---- figure 1: map and across-fault section of each catalog (same events, same scales)
    W = 5.2                                             # panel width (cm)
    reg = [-121.27, -121.12, 36.545, 36.650]
    xr, zr = [-4.0, 2.5], [2.0, 10.5]                   # across-fault section, true scale
    Hs = W * (zr[1] - zr[0]) / (xr[1] - xr[0])
    fig = pygmt.Figure()
    pygmt.config(FONT="8p,Helvetica", MAP_FRAME_TYPE="plain", FORMAT_GEO_MAP="ddd.xx")
    gap = 1.25
    for j, k in enumerate(("ncsn", "ddrt", "ours")):
        P = X[k]
        lon, lat = ll(P)
        if j:
            fig.shift_origin(xshift=f"{W + gap}c", yshift=f"{Hs + 1.5}c")
        fig.basemap(region=reg, projection=f"M{W}c", frame=["WSne" if j == 0 else "wSne", "xa0.05", "ya0.05"])
        if j == 0:
            fig.basemap(map_scale="jBL+w2k+o0.25c/0.3c+c36.6")
        for ln in segs:
            fig.plot(x=ln[:, 0], y=ln[:, 1], pen="1.1p,brown")
        fig.plot(x=lon, y=lat, style="c", size=size(mag), fill=COL[k], pen="0.2p,black")
        fig.plot(x=[lon[h]], y=[lat[h]], style="a0.36c", fill=RED, pen="0.4p,black")
        title(fig, NAME[k])
        letter(fig, "(" + "abc"[j] + ")", -1.15 if j == 0 else -0.35)
        # across-fault section (events along the section stretch of the line)
        fig.shift_origin(yshift=f"-{Hs + 1.5}c")
        a = (P[:, :2] - P[h, :2]) @ fr.along_u
        sec = (a >= g["section"][0]) & (a <= g["section"][1])
        ne = m[f"{k}_ne_of_trace_km"].values
        fig.basemap(region=xr + zr, projection=f"X{W}c/-{Hs:.3f}c",
                    frame=["WSne" if j == 0 else "wSne", "xa1f0.5+lDistance NE of the fault trace (km)",
                           "ya2f1+lDepth (km)" if j == 0 else "ya2f1"])
        fig.plot(x=[0, 0], y=zr, pen="0.8p,brown,4_2")
        # planes through the catalog's own mainshock: the moment tensors (dashed) and the aftershocks' fit (solid)
        z = np.array(zr)
        for (lab, (s, dp)), pen in zip(MT_PLANES.items(), ("0.9p,gray30,6_3", "1.1p,gray30,1_2")):
            off = np.radians(s - STRIKE)
            tan_app = np.tan(np.radians(dp)) * np.cos(off)          # apparent dip in the section
            fig.plot(x=ne[h] - (z - P[h, 2]) / tan_app, y=z, pen=pen)   # dips SW (right-hand rule, strike ~137)
        st, dp, n, c = plane_fit(P[aft])
        n_ne, n_al = n[:2] @ fr.ne_u, n[:2] @ fr.along_u
        cx = float(np.mean(ne[aft]))
        fig.plot(x=cx - n[2] / n_ne * (z - c[2]), y=z, pen="1.1p,black")
        fig.plot(x=ne[sec], y=P[sec, 2], style="c", size=size(mag[sec]), fill=COL[k], pen="0.2p,black")
        fig.plot(x=[ne[h]], y=[P[h, 2]], style="a0.36c", fill=RED, pen="0.4p,black")
        r = res["catalogs"][k]
        box = dict(justify="LB", font="7p,Helvetica", fill="white", clearance="0.05c/0.04c")
        fig.text(x=xr[0] + 0.12, y=zr[1] - 0.15, **box,
                 text=f"Fitted dip {r['plane_dip']:.0f}° ({r['plane_dip_95pct'][0]:.0f}-{r['plane_dip_95pct'][1]:.0f})")
        fig.text(x=xr[0] + 0.12, y=zr[1] - 0.7, **box,
                 text=f"Width {r['width_m']} m; {r['off_line_gt_0p5km']} > 0.5 km off")
        fig.text(x=xr[0] + 0.15, y=zr[0] + 0.2, justify="LT", font="7p,Helvetica,gray30", text="SW")
        fig.text(x=xr[1] - 0.15, y=zr[0] + 0.2, justify="RT", font="7p,Helvetica,gray30", text="NE")
        letter(fig, "(" + "def"[j] + ")", -1.15 if j == 0 else -0.35)
    key = Path(tmp) / "key1.txt"
    key.write_text("\n".join([
        "N 3",
        f"S 0.2c c 0.15c {COL['ncsn']} 0.2p,black 0.4c NCSN catalog",
        f"S 0.2c c 0.15c {COL['ddrt']} 0.2p,black 0.4c LDEO real-time DD",
        f"S 0.2c c 0.15c {COL['ours']} 0.2p,black 0.4c PocketQuake-global",
        f"S 0.2c a 0.3c {RED} 0.4p,black 0.4c Mw 4.98 mainshock",
        "S 0.2c - 0.4c - 1.1p,brown 0.45c San Andreas fault (NSHM23)",
        "S 0.2c - 0.4c - 1.1p,black 0.45c Plane fitted to the aftershocks",
        "S 0.2c - 0.4c - 0.9p,gray30,6_3 0.45c NC moment tensor (137/73)",
        "S 0.2c - 0.4c - 1.1p,gray30,1_2 0.45c USGS Mww (135/82)"]) + "\n")
    fig.shift_origin(xshift=f"-{2 * (W + gap)}c", yshift=f"{Hs + 1.5}c")
    fig.basemap(region=reg, projection=f"M{W}c", frame="+n")             # the legend sits on map (a)
    fig.legend(spec=str(key), position=f"JTL+jBL+o0c/0.75c+w{3 * W + 2 * gap:.2f}c", box=False)
    out1 = Path(out_dir) / "ldeo_ddrt_comparison.png"
    fig.savefig(out1, dpi=250)
    fig.savefig(out1.with_suffix(".pdf"))

    # ---- figure 2: ours and DDRT after removing the common shift, tied event by event (strike frame, true scale)
    Pd, Po = X["ddrt"] - g["shift"], X["ours"]
    A = lambda P: np.c_[(P[:, :2] - Po[h, :2]) @ fr.along_u, (P[:, :2] - Po[h, :2]) @ fr.ne_u, P[:, 2]]
    ad, ao = A(Pd), A(Po)
    ar, cr, zr2 = [-4.0, 7.0], [-1.5, 1.5], [3.0, 9.5]
    W2 = 15.0
    fig = pygmt.Figure()
    pygmt.config(FONT="8p,Helvetica", MAP_FRAME_TYPE="plain")
    H1 = W2 * (cr[1] - cr[0]) / (ar[1] - ar[0])
    H2 = W2 * (zr2[1] - zr2[0]) / (ar[1] - ar[0])
    for i, (yi, yr, lab, flip) in enumerate(((1, cr, "Across strike, NE + (km)", False), (2, zr2, "Depth (km)", True))):
        if i:
            fig.shift_origin(yshift=f"-{H2 + 1.6}c")
        fig.basemap(region=ar + yr, projection=f"X{W2}c/{'-' if flip else ''}{H1 if i == 0 else H2:.3f}c",
                    frame=["WSne", "xa1f0.5+lAlong strike 137° from the mainshock (km)" if i else "xa1f0.5",
                           f"ya1f0.5+l{lab}"])
        for p, o in zip(ad, ao):
            fig.plot(x=[o[0], p[0]], y=[o[yi], p[yi]], pen="0.4p,gray55")
        fig.plot(x=ad[:, 0], y=ad[:, yi], style="c", size=size(mag), fill=COL["ddrt"], pen="0.2p,black")
        fig.plot(x=ao[:, 0], y=ao[:, yi], style="c", size=size(mag), fill=COL["ours"], pen="0.2p,black")
        fig.plot(x=[0], y=[ao[h, yi]], style="a0.36c", fill=RED, pen="0.4p,black")
        letter(fig, "(" + "ab"[i] + ")", -1.15)
    rel = res["relative_after_shift"]
    sh = res["shift_ddrt_minus_ours_m"]
    key = Path(tmp) / "key2.txt"
    key.write_text("\n".join([
        "N 3.2 4.8 2.4 3.2",
        f"S 0.2c c 0.15c {COL['ours']} 0.2p,black 0.4c PocketQuake-global",
        f"S 0.2c c 0.15c {COL['ddrt']} 0.2p,black 0.4c LDEO real-time DD (shifted)",
        "S 0.2c - 0.4c - 0.6p,gray55 0.45c Same event",
        f"S 0.2c a 0.3c {RED} 0.4p,black 0.4c Mainshock (ours)"]) + "\n")
    fig.shift_origin(yshift=f"{H2 + 1.6}c")
    fig.basemap(region=ar + cr, projection=f"X{W2}c/{H1:.3f}c", frame="+n")       # the legend sits on panel (a)
    fig.legend(spec=str(key), position=f"JTL+jBL+o0c/0.75c+w{W2}c", box=False)
    fig.text(text=(f"LDEO positions shifted by {-sh['east']:+d} m east, {-sh['north']:+d} m north, {-sh['down']:+d} m down; "
                   f"then median |difference| {rel['along_strike']['median_m']} m along, "
                   f"{rel['across_strike']['median_m']} m across, {rel['depth']['median_m']} m in depth"),
             position="TL", justify="BL", offset="0c/0.2c", no_clip=True, font="7p,Helvetica,gray30")
    out2 = Path(out_dir) / "ldeo_ddrt_differences.png"
    fig.savefig(out2, dpi=250)
    fig.savefig(out2.with_suffix(".pdf"))
    return out1, out2


def main():
    kcr = HERE.parents[1] / "external" / "PocketQuake" / "external" / "korea-cluster-relocation"
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reloc", default=str(kcr / "pipeline" / "runs" / "pinnacles2026" / "2.HypoDD" / "02.dt.cc" / "hypoDD.reloc"),
                    help="hypoDD.reloc, or a CSV with event, lat, lon, depth (and solution)")
    ap.add_argument("--solution", default="hypodd_cc", help="with a CSV that has a `solution` column")
    ap.add_argument("--catalog", default=str(HERE / "comcat_sequence.csv"))
    ap.add_argument("--ddrt", default=str(HERE / "ldeo_ddrt_sequence.txt"))
    ap.add_argument("--trace", default=str(HERE / "usgs_nshm23_san_andreas.geojson"))
    ap.add_argument("--classes", default=None, help="CSV with id, class (aftershock = 'aftershock'); default: "
                    f"after the mainshock and within {ZONE_KM:g} km of it")
    ap.add_argument("--section", default="-4,7", help="along-strike stretch (km) for the section and the widths")
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1128)
    ap.add_argument("--out-dir", default=str(HERE.parents[1] / "work" / "pinnacles2026" / "ldeo_ddrt"))
    ap.add_argument("--no-figures", action="store_true")
    a = ap.parse_args()
    cat = pd.read_csv(a.catalog)
    ours = read_ours(a.reloc, cat, a.solution)
    cls = pd.read_csv(a.classes).set_index("id")["class"] if a.classes else None
    segs = trace_segments(a.trace)
    res, m, g = compare(ours, cat, read_ddrt(a.ddrt), segs, cls, tuple(float(x) for x in a.section.split(",")),
                        a.bootstrap, a.seed)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    res["inputs"] = dict(reloc=str(a.reloc), ddrt=str(a.ddrt), catalog=str(a.catalog), classes=a.classes)
    json.dump(res, open(out / "ldeo_ddrt_comparison.json", "w"), indent=1)
    m.round(5).to_csv(out / "ldeo_ddrt_comparison.csv", index=False)
    print(json.dumps({k: v for k, v in res.items() if k != "inputs"}, indent=1))
    if not a.no_figures:
        for f in figures(res, m, g, segs, out):
            print(f)


if __name__ == "__main__":
    main()
