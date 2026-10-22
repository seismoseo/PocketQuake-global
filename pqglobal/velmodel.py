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
