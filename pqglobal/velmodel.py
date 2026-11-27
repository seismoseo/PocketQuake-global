"""Gradient velocity models -> the constant-velocity layers HYPOINVERSE, HypoDD and SKHASH take.

HYPOINVERSE's layered (CRH) models hold at most 20 layers including the half-space. A gradient
model (velocity linear between depth nodes, e.g. NCSN CON) is cut into at most that many layers:
layers of equal thickness inside each gradient segment, each with the velocity that keeps its
vertical travel time (thickness / integral of dz/v), and the number of layers per segment chosen
to minimise the largest first-arrival P error over the distances and depths that matter.

The error is measured against the gradient model itself, sampled in 10 m layers, with the same
flat-earth first-arrival solver (direct ray or head wave along any interface below the source, the
earliest one; the station on the model surface, as HYPOINVERSE puts it). With one Vp/Vs ratio
the S error is the P error times Vp/Vs.
"""
from __future__ import annotations

import itertools

import numpy as np

# NCSN "CON" (Concord-Calaveras): depth below the model surface (km) and Vp (km/s), linear in
# between, constant below the last node. Oppenheimer et al. (1993) OFR 93-578; NCEDC
# shadow2000.pdf Table 8; NCSN cor.crt (2023).
CON_NODES = ((0.0, 2.50), (2.0, 4.70), (12.0, 5.77), (24.0, 6.12), (26.0, 7.95))
MAX_LAYERS = 20           # HYPOINVERSE CRH limit, half-space included
_trapezoid = getattr(np, "trapezoid", None) or np.trapz     # numpy >= 2 renamed trapz


def _v_at(nodes, z):
    zs, vs = zip(*nodes)
    return np.interp(z, zs, vs)


def fine_model(nodes, dz=0.01, zmax=40.0):
    """The gradient model as dz-thick layers (velocity at each layer's middle) + half-space."""
    tops = np.arange(0.0, nodes[-1][0], dz)
    v = _v_at(nodes, tops + dz / 2)
    return np.append(tops, nodes[-1][0]), np.append(v, nodes[-1][1])


def cut(nodes, counts):
    """Layer tops and velocities: counts[i] equal layers in gradient segment i, + half-space.
    Each layer keeps its vertical travel time (slowness averaged over the layer)."""
    tops, vel = [], []
    for (z0, v0), (z1, v1), n in zip(nodes[:-1], nodes[1:], counts):
        edges = np.linspace(z0, z1, n + 1)
        for a, b in zip(edges[:-1], edges[1:]):
            zz = np.linspace(a, b, 201)
            slow = _trapezoid(1.0 / _v_at(nodes, zz), zz) / (b - a)
            tops.append(a)
            vel.append(1.0 / slow)
    tops.append(nodes[-1][0])
    vel.append(nodes[-1][1])
    return np.array(tops), np.array(vel)


def first_arrival(tops, v, zs, X):
    """Earliest P time (s) from a source at depth zs (km) to surface stations at offsets X (km)."""
    X = np.asarray(X, float)
    bottoms = np.append(tops[1:], np.inf)
    above = np.clip(np.minimum(bottoms, zs) - tops, 0, None)          # thickness above source
    m = above > 0
    th_a, v_a = above[m], v[m]
    # direct (up-going) ray, bisection on the ray parameter
    pmax = 1.0 / v_a.max()
    lo, hi = np.zeros_like(X), np.full_like(X, pmax * (1 - 1e-12))
    for _ in range(60):
        p = (lo + hi) / 2
        pv = p[:, None] * v_a[None, :]
        x = (th_a * pv / np.sqrt(1 - pv ** 2)).sum(1)
        lo, hi = np.where(x < X, p, lo), np.where(x < X, hi, p)
    pv = p[:, None] * v_a[None, :]
    t_best = (th_a / (v_a * np.sqrt(1 - pv ** 2))).sum(1)
    t_best = np.where(X == 0, (th_a / v_a).sum(), t_best)
    # head waves along the top of every layer at or below the source, faster than all above
    for k in range(1, len(tops)):
        if tops[k] < zs or v[k] <= v[:k].max():
            continue
        pk = 1.0 / v[k]
        up = np.clip(np.minimum(bottoms[:k], tops[k]) - tops[:k], 0, None)       # surface..z_k
        down = np.clip(np.minimum(bottoms[:k], tops[k]) - np.maximum(tops[:k], zs), 0, None)
        th = up + down
        eta = np.sqrt(1.0 / v[:k] ** 2 - pk ** 2)
        xc = (th * pk * v[:k] / np.sqrt(1 - (pk * v[:k]) ** 2)).sum()
        t = X * pk + (th * eta).sum()
        t_best = np.where(X >= xc, np.minimum(t_best, t), t_best)
    return t_best


def reference_times(nodes, depths, dists):
    """First-arrival P times (s) of the gradient model, array [depth, dist]."""
    ft, fv = fine_model(nodes)
    return np.array([first_arrival(ft, fv, z, dists) for z in depths])


def errors(tops, vel, depths, dists, ref):
    """P time of the layered model minus the gradient model's (`ref`), ms, array [depth, dist]."""
    return np.array([first_arrival(tops, vel, z, dists) for z in depths]) * 1e3 - ref * 1e3


def best_cut(nodes, depths, dists, max_layers=MAX_LAYERS):
    """(max |error| ms, per-segment layer counts summing to max_layers - 1) with the smallest
    largest P error over depths x dists."""
    ref = reference_times(nodes, depths, dists)
    best = None
    for counts in itertools.product(range(1, max_layers), repeat=len(nodes) - 1):
        if sum(counts) != max_layers - 1:
            continue
        tops, vel = cut(nodes, counts)
        e = np.abs(errors(tops, vel, depths, dists, ref)).max()
        if best is None or e < best[0]:
            best = (e, counts)
    return best


def greedy_cut(nodes, depths, dists, max_layers=MAX_LAYERS):
    """Like best_cut, but greedy (for models with many gradient segments): one layer per segment,
    then one layer at a time to the segment whose split lowers the largest P error most."""
    ref = reference_times(nodes, depths, dists)
    counts = [1] * (len(nodes) - 1)

    def err(c):
        tops, vel = cut(nodes, c)
        return np.abs(errors(tops, vel, depths, dists, ref)).max()
    while sum(counts) < max_layers - 1:
        trials = []
        for i in range(len(counts)):
            c = list(counts)
            c[i] += 1
            trials.append((err(c), c))
        counts = min(trials, key=lambda t: t[0])[1]
    return err(counts), tuple(counts)


# ----------------------------------------------------------------------- NCSN multiple models
def read_crt_nodes(path):
    """((depth_km, vp_km_s), ...) of an NCSN gradient model file (.crt, as read by HYPOINVERSE's CRT
    command): header line 2 holds the node depths, line 3 the velocities, F5.2 each."""
    with open(path) as f:
        lines = f.read().splitlines()
    fields = lambda s: [float(s[i:i + 5]) for i in range(0, len(s.rstrip()), 5) if s[i:i + 5].strip()]
    z, v = fields(lines[1]), fields(lines[2])
    if len(z) != len(v) or z[0] != 0.0:
        raise ValueError(f"{path}: cannot read the node lines of this .crt file")
    return tuple(zip(z, v))


def ncsn_pair(multmod_path, code):
    """NCSN's setup of model `code` (e.g. 'gab') and its alternate, from its multiple-model command
    file (multmod5a.hyp): their .crt files, delay files, P-to-S ratios and the alternate-station
    file, exactly as NCSN reads them. Paths are resolved against the file's folder."""
    import os
    import re
    base = os.path.dirname(os.path.abspath(multmod_path))
    crt, dly, psm, alt, altfile = {}, {}, {}, {}, None
    for raw in open(multmod_path):
        m = re.match(r"\s*(crt|del|psm|alt)\s+(-?\d+)\s+(.*)$", raw, re.IGNORECASE)
        if not m:                                     # comments (*), nod, mul, pos, ...
            continue
        cmd, n, rest = m.group(1).lower(), int(m.group(2)), m.group(3)
        if cmd in ("crt", "del"):
            path = os.path.normpath(os.path.join(base, re.search(r"'([^']+)'", rest).group(1)))
            if cmd == "crt":
                crt[n] = path
            elif n == -1:
                altfile = path
            else:
                dly[n] = path
        elif cmd == "psm":
            x = rest.split()
            psm[n] = (float(x[0]), float(x[1]))
        else:
            alt[n] = int(rest.split()[0])
    prim = [n for n, p in crt.items() if os.path.basename(p).lower() == f"{code.lower()}.crt"]
    if not prim or prim[0] not in alt:
        raise ValueError(f"{multmod_path}: no model '{code}' with an alternate (ALT) found")
    p, a = prim[0], alt[prim[0]]
    return dict(primary=dict(crt=crt[p], delays=dly.get(p), psm=psm.get(p, (1.75, 0.0))),
                alternate=dict(crt=crt[a], delays=dly.get(a), psm=psm.get(a, (1.75, 0.0))),
                alternate_stations=altfile)


def slowness_average(nodes_a, nodes_b):
    """Gradient-model nodes at the union of both models' node depths, velocity = harmonic mean
    (the slowness average) of the two models there; constant below the last node."""
    z = sorted({d for d, _ in nodes_a} | {d for d, _ in nodes_b})
    va, vb = _v_at(nodes_a, np.array(z)), _v_at(nodes_b, np.array(z))
    return tuple((float(d), float(round(2 / (1 / x + 1 / y), 4))) for d, x, y in zip(z, va, vb))
