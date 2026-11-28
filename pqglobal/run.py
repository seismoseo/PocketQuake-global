"""pocketquake-global: from a USGS catalog to PocketQuake's relocation chain, via FDSN.

Stages (each one writes into --work and is reused on a rerun):
  1. catalog       USGS ComCat query (box or circle, UTC) or a ComCat CSV -> event_catalog.csv
                   (--largest N keeps the N largest events, e.g. for a validation subset)
  2. picks         the network's own picks for those events (FDSN event service) -> network_picks.csv
  3. stations      the stations the network picked, within --radius-km, with a three-component
                   sensor at >= --min-sps -> station_table/ (FDSN_station.csv, sensors.csv,
                   inventory.xml)
  4. run           PocketQuake's orchestrator: scaffold the cluster, fetch the event waveforms
                   (FDSN dataselect), picking, HYPOINVERSE, HypoDD dt.ct/dt.cc, SKHASH, notebook;
                   pairs with an event of magnitude >= --onset-mag (2.0) are cross-correlated on the
                   onset only

--stop-after stations shows the station set before anything is downloaded. Defaults: dry run
(prints what it would do); --apply acts.
"""
from __future__ import annotations

import argparse
import math
import os
import re
import shutil

import pandas as pd

import pqglobal                                    # noqa: F401  (sets the PocketQuake paths)
from pqglobal import catalog as catmod
from pqglobal import network_picks, stations, velmodel
from pqglobal.waveforms import FDSNFetcher, event_id

STATION_FILES = ("FDSN_station.csv", "sensors.csv", "inventory.xml")


def _rotation_reference(value):
    """Parse the optional shared source point used for R/T back azimuths."""
    try:
        lat, lon = (float(x) for x in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected LAT,LON") from exc
    if not (math.isfinite(lat) and math.isfinite(lon) and
            -90 <= lat <= 90 and -180 <= lon <= 180):
        raise argparse.ArgumentTypeError("LAT,LON must be finite and within geographic bounds")
    return lat, lon


def _velocity_model_arg(name, counts, vpvs, layers_csv=None):
    """Source text of the cluster module's velocity model: the NCSN CON gradient model cut into
    layers (`counts` per gradient segment), or the layers of `layers_csv` (top_km, vp_km_s)."""
    if layers_csv:
        lay = pd.read_csv(layers_csv)
        if not {"top_km", "vp_km_s"} <= set(lay.columns):
            raise SystemExit(f"{layers_csv}: need columns top_km, vp_km_s")
        if len(lay) > velmodel.MAX_LAYERS:
            raise SystemExit(f"{layers_csv}: {len(lay)} layers; HYPOINVERSE takes at most "
                             f"{velmodel.MAX_LAYERS} (half-space included)")
        tops, vp = lay["top_km"].to_numpy(float), lay["vp_km_s"].to_numpy(float)
    else:
        tops, vp = velmodel.cut(velmodel.CON_NODES, counts)
    t = ", ".join(f"{x:.3f}" for x in tops)
    v = ", ".join(f"{x:.3f}" for x in vp)
    return f'layered_velmodel("{name}", ({t}), ({v}), {vpvs})'


def _ncsn_velocity_model_arg(name, code, multmod, alt_extra, work, vpvs=None):
    """Source text of the cluster module's velocity model following NCSN's own setup for a region
    with two models across a fault (e.g. 'gab': Gabilan + its alternate Diablo, Bear Valley):
    HYPOINVERSE reads NCSN's two gradient models (.crt), their station delays and P-to-S ratios, and
    the alternate-station list (NCSN's file + `alt_extra`, lines 'STA   NET  A') exactly as NCSN
    does (CRT/DEL/PSM/ALT, DEL -1); HypoDD gets the slowness average of the two cut into <= 20 layers;
    SKHASH traces takeoffs through that average and both NCSN models. Returns (text, vpvs, note)."""
    import numpy as np
    pair = velmodel.ncsn_pair(multmod, code)
    prim, alt = pair["primary"], pair["alternate"]
    gn, an = velmodel.read_crt_nodes(prim["crt"]), velmodel.read_crt_nodes(alt["crt"])
    avg = velmodel.slowness_average(gn, an)
    err, counts = velmodel.greedy_cut(avg, np.arange(1, 16, 1.0), np.arange(0, 121, 5.0))
    tops, vp = velmodel.cut(avg, counts)
    altf = os.path.join(os.path.abspath(work), "ncsn_alternate_stations.del")
    with open(altf, "w") as f:
        f.write(open(pair["alternate_stations"]).read().rstrip("\n") + "\n")
        if alt_extra:
            f.write(open(alt_extra).read().rstrip("\n") + "\n")
    b = os.path.basename
    lines = (f"DEL -1 '{{dir}}/{b(altf)}'",
             f"CRT 1 '{{dir}}/{b(prim['crt'])}'", f"DEL 1 '{{dir}}/{b(prim['delays'])}'",
             f"CRT 2 '{{dir}}/{b(alt['crt'])}'", f"DEL 2 '{{dir}}/{b(alt['delays'])}'",
             f"PSM 1 {prim['psm'][0]} {prim['psm'][1]}", f"PSM 2 {alt['psm'][0]} {alt['psm'][1]}",
             "ALT 1 2")
    files = tuple(dict.fromkeys(os.path.abspath(x) for x in
                                (prim["crt"], alt["crt"], prim["delays"], alt["delays"], altf)))
    vpvs = vpvs or prim["psm"][0]
    t = ", ".join(f"{x:.3f}" for x in tops)
    v = ", ".join(f"{x:.3f}" for x in vp)
    tk = ", ".join(f"(({', '.join(str(z) for z, _ in nd)}), ({', '.join(str(w) for _, w in nd)}))"
                   for nd in (gn, an))
    text = (f'layered_velmodel("{name}", ({t}), ({v}), {vpvs}, hyp_lines={lines!r}, '
            f'hyp_files={files!r}, takeoff_models=({tk},))')
    note = (f"NCSN {code.upper()} + alternate {b(alt['crt'])}; HypoDD/SKHASH average cut {counts} "
            f"(largest P error {err:.1f} ms, 1-15 km deep, 0-120 km)")
    return text, vpvs, note


ONSET_WINDOW = dict(pre=0.05, post=0.05, bandpass=(1, 40))


def onset_events(cat, onset_mag):
    """Event ids (the pipeline's UTC second) of the catalog events with magnitude >= onset_mag, largest
    first; () for 'none'. Pairs that include one of them are cross-correlated on the onset only."""
    if onset_mag is None or str(onset_mag).lower() == "none":
        return ()
    big = cat[cat.magnitude >= float(onset_mag)].sort_values("magnitude", ascending=False, kind="mergesort")
    return tuple(event_id(t) for t in big.time)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--name", required=True, help="cluster name, e.g. sanramon2015")
    ap.add_argument("--region", required=True, help="display name, e.g. SanRamon2015")
    ap.add_argument("--work", required=True, help="working directory for stages 1-3")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--usgs-box", help="lat0,lat1,lon0,lon1 (ComCat query)")
    g.add_argument("--catalog-csv", help="a ComCat CSV export")
    ap.add_argument("--start"), ap.add_argument("--end")
    ap.add_argument("--minmag", type=float)
    ap.add_argument("--largest", type=int, help="keep only the N largest events")
    ap.add_argument("--picks-service", help="FDSN event service with picks (default from ids)")
    ap.add_argument("--data-service", default="NCEDC", help="FDSN station/dataselect service")
    ap.add_argument("--radius-km", type=float, default=60.0)
    ap.add_argument("--min-sps", type=float, default=100.0)
    ap.add_argument("--velmodel", required=True,
                    help="'con' (NCSN Concord-Calaveras, Northern California; cut to 20 layers), a "
                         "CSV of layers with columns top_km, vp_km_s (<= 20 rows, half-space last), or "
                         "'ncsn:<code>' for an NCSN model with an alternate across a fault, e.g. "
                         "'ncsn:gab' (Bear Valley), read from --ncsn-multmod")
    ap.add_argument("--ncsn-multmod", help="NCSN multiple-model file (multmod5a.hyp of NCEDC's "
                                           "hypoinv-standalone package), for --velmodel ncsn:<code>")
    ap.add_argument("--alternate-extra", help="extra alternate-model stations ('STA   NET  A' lines) "
                                              "for --velmodel ncsn:<code>, e.g. stations newer than NCSN's list")
    ap.add_argument("--vpvs", type=float, default=None,
                    help="one Vp/Vs for the S model and HypoDD (default 1.747; NCSN's PSM value with ncsn:)")
    ap.add_argument("--con-layers", default="5,8,5,1",
                    help="CON layers per gradient segment (+ half-space; 20 in all at most)")
    ap.add_argument("--pre", type=float, default=240.0)
    ap.add_argument("--post", type=float, default=360.0)
    ap.add_argument("--fetch-workers", type=int, default=2)
    ap.add_argument("--cores", type=int, default=2, help="xcorr worker cap (memory!)")
    ap.add_argument("--onset-mag", default="2.0",
                    help="cross-correlate every event pair that includes an event of at least this "
                         "magnitude on the onset only (+-0.05 s around the pick, 1-40 Hz, instead of "
                         "+-0.5 s, 5-20 Hz), so the larger events are placed by where their rupture "
                         "started; 'none' turns it off (default 2.0)")
    ap.add_argument("--s-on-t", action="store_true",
                    help="measure S differential times on T after N/E to R/T rotation; P stays on Z")
    ap.add_argument("--orientation-csv", metavar="PATH",
                    help="time-dependent residual N/E misorientation table (implies --s-on-t); "
                         "FDSN metadata orientation is already applied")
    ap.add_argument("--orientation-missing", choices=("error", "zero", "skip"), default=None,
                    help="policy for missing orientation epochs (default: error)")
    ap.add_argument("--rotation-reference", type=_rotation_reference, metavar="LAT,LON",
                    help="shared source point for R/T back azimuths (default: each event)")
    ap.add_argument("--stop-after", choices=("catalog", "picks", "stations", "run"), default="run")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args(argv)
    if (not a.s_on_t and not a.orientation_csv
            and (a.orientation_missing is not None or a.rotation_reference is not None)):
        ap.error("orientation settings require --s-on-t or --orientation-csv")
    os.makedirs(a.work, exist_ok=True)

    # 1. catalog
    cat_csv = os.path.join(a.work, "event_catalog.csv")
    if not os.path.exists(cat_csv):
        if a.usgs_box:
            box = tuple(float(x) for x in a.usgs_box.split(","))
            df = catmod.from_usgs(a.start, a.end, box=box, minmagnitude=a.minmag,
                                  cache_csv=os.path.join(a.work, "comcat_raw.csv"))
        else:
            df = catmod.from_usgs_csv(a.catalog_csv)
        if a.largest:
            df = df.sort_values("magnitude", ascending=False, kind="mergesort").head(a.largest)
        print("[catalog]", catmod.write(df, cat_csv))
    cat = pd.read_csv(cat_csv)
    (lat, lon), bounds = catmod.centroid_and_bounds(cat)
    print(f"[catalog] {len(cat)} events, M {cat.magnitude.min():.2f}-{cat.magnitude.max():.2f}, "
          f"centroid {lat}, {lon}")
    if a.stop_after == "catalog":
        return

    # 2. network picks
    picks_csv = os.path.join(a.work, "network_picks.csv")
    if not os.path.exists(picks_csv):
        print("[picks]", network_picks.fetch(cat, picks_csv, os.path.join(a.work, "picks_cache"),
                                             service=a.picks_service))
    if a.stop_after == "picks":
        return

    # 3. stations
    st_dir = os.path.join(a.work, "station_table")
    if not all(os.path.exists(os.path.join(st_dir, f)) for f in STATION_FILES):
        t = pd.to_datetime(cat["time"])
        print("[stations]", stations.select_picked(
            picks_csv, (lat, lon), str(t.min().floor("D").isoformat()),
            str((t.max().floor("D") + pd.Timedelta(days=1)).isoformat()), st_dir,
            service=a.data_service, radius_km=a.radius_km, min_sps=a.min_sps))
    sens = pd.read_csv(os.path.join(st_dir, "sensors.csv"), dtype={"location": str})
    used = sens[sens.verdict == "used"]
    bands = sorted(used.band.unique(), key=stations._rank)
    print(f"[stations] {len(used)} stations; sensors by band {used.band.value_counts().to_dict()}")
    if a.stop_after == "stations":
        return

    # 4. PocketQuake orchestrator
    big = cat.loc[cat.magnitude.idxmax()]
    if a.velmodel.lower().startswith("ncsn:"):
        if not a.ncsn_multmod:
            raise SystemExit("--velmodel ncsn:<code> needs --ncsn-multmod")
        code = a.velmodel.split(":", 1)[1]
        vm_name = f"ncsn_{code.lower()}"[:16]
        vm_arg, a.vpvs, note = _ncsn_velocity_model_arg(vm_name, code, a.ncsn_multmod, a.alternate_extra,
                                                        a.work, a.vpvs)
        print("[velmodel]", note)
    else:
        a.vpvs = a.vpvs or 1.747
        layers_csv = None if a.velmodel.lower() == "con" else a.velmodel
        vm_name = ("con" if layers_csv is None else       # names HYPOINVERSE's <name>_p.crh files
                   re.sub(r"[^A-Za-z0-9_]", "_", os.path.splitext(os.path.basename(layers_csv))[0])[:16])
        counts = tuple(int(x) for x in a.con_layers.split(","))
        vm_arg = _velocity_model_arg(vm_name, counts, a.vpvs, layers_csv)
    # the station table already is the selection (stage 3); the cluster's own radius cut, measured
    # from this catalog's centroid (a subset's differs), gets a margin so it never cuts again
    factory_args = (f"velocity_models=({vm_arg},), "
                    f"sensors={tuple(bands)!r}, radius_km={a.radius_km + 10}, hypodd_vpvs={a.vpvs}, "
                    f'mainshock_event_id="{event_id(big.time)}", '
                    f'magnitude_label="USGS ComCat magnitude"')
    onset = onset_events(cat, a.onset_mag)
    if onset:
        factory_args += (f", xcorr_pair_overrides={{frozenset({onset!r}): "
                         f"dict(pre={ONSET_WINDOW['pre']}, post={ONSET_WINDOW['post']}, "
                         f"bandpass={ONSET_WINDOW['bandpass']!r})}}")

    def station_tables_fn(dst):
        os.makedirs(dst, exist_ok=True)
        for f in STATION_FILES:
            shutil.copyfile(os.path.join(st_dir, f), os.path.join(dst, f))

    fetcher = FDSNFetcher(service=a.data_service, pre_s=a.pre, post_s=a.post,
                          workers=a.fetch_workers)
    kw = dict(region=a.region, wf_backend="external", fetcher=fetcher, velmodel=vm_name,
              cores=a.cores, s_on_t=a.s_on_t or bool(a.orientation_csv),
              orientation_csv=a.orientation_csv, orientation_missing=a.orientation_missing or "error",
              rotation_reference=a.rotation_reference,
              spec_extra=dict(factory="global_cluster",
                                             factory_imports=("layered_velmodel",),
                                             factory_args=factory_args,
                                             station_tables_fn=station_tables_fn,
                                             station_master_files=("FDSN_station.csv",)))
    print(f"[run] cluster {a.name} ({a.region}) epicenter {lat}, {lon}; bounds {bounds}\n"
          f"      {factory_args}")
    print(f"[run] onset window for {len(onset)} events of magnitude >= {a.onset_mag}" if onset
          else "[run] no onset window (--onset-mag none, or no event that large)")
    if kw["s_on_t"]:
        print(f"[run] S on T; residual orientation table={a.orientation_csv or 'none (0 deg)'}; "
              f"missing={kw['orientation_missing']}; reference={a.rotation_reference or 'each event'}")
    if not a.apply:
        print("[run] dry run: pass --apply to scaffold, download and run the pipeline")
        return
    from pocketquake.orchestrate import orchestrate
    print(orchestrate(cat_csv, a.name, (lat, lon), bounds, **kw))


if __name__ == "__main__":
    main()
