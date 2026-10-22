"""Event waveforms from an FDSN dataselect service into the pipeline's per-event SAC layout.

For every catalog event (folder name = the origin's UTC second, YYYYMMDDHHMMSS, as the pipeline
names events) and every sensor in station_table/inventory.xml:

  1. fetch origin - pre_s .. origin + post_s (default -240 s .. +360 s, as the KMA/KIGAM event
     windows), bulk requests of BULK_CHUNK sensors, at most `workers` events at a time;
  2. per sensor, keep one gap-free stretch shared by all three components (the one that covers
     the origin, else the longest);
  3. divide each channel by its signed overall sensitivity (counts -> m/s or m/s^2), so a
     channel wired with negative gain gets the right polarity;
  4. rotate to Z/N/E with the metadata's azimuth and dip: NCSN/BK borehole channels named 1/2/3
     (1 vertical), misoriented horizontals and downward-pointing verticals all come out as
     ground motion up / north / east;
  5. write SAC (float32, sampling rate unchanged; the pipeline resamples where it needs to) with
     station and event coordinates, lcalda=1, reference time = origin, into
     <event>/FDSN.{v,a}/SAC/<band>/NET.STA.LOC.CHA.D.YYYY.DDD.HHMMSS.SAC
     (.a = accelerometer, instrument code N or G; .v = everything else).

Each event is built in <waveforms_dir>/.staging/ and moved in place when complete, then marked
by <event>/FDSN.done.json (counts and the sensors without data), so a rerun skips finished
events and redoes interrupted ones.
"""
from __future__ import annotations

import json
import os
import shutil
import time as _time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
from obspy import Stream, UTCDateTime, read_inventory

BULK_CHUNK = 40
MIN_COMMON_S = 30.0          # a sensor whose three components share less than this is skipped
DONE = "FDSN.done.json"


def event_id(t) -> str:
    return UTCDateTime(t).strftime("%Y%m%d%H%M%S")


def _sensors(inv):
    """[(net, sta, loc, band)] of the inventory, one entry per sensor."""
    out = set()
    for n in inv:
        for s in n:
            for c in s:
                out.add((n.code, s.code, c.location_code, c.code[:2]))
    return sorted(out)


def _common_stretch(st3, origin):
    """One gap-free window shared by the components: per component the piece covering the origin
    (else the longest), then their overlap. Returns a trimmed Stream or None."""
    pieces = {}
    for tr in st3.split():
        pieces.setdefault(tr.stats.channel, []).append(tr)
    if len(pieces) < 3:
        return None
    pick = []
    for ch, trs in pieces.items():
        cover = [t for t in trs if t.stats.starttime <= origin <= t.stats.endtime]
        pick.append(max(cover or trs, key=lambda t: t.stats.npts))
    t0 = max(t.stats.starttime for t in pick)
    t1 = min(t.stats.endtime for t in pick)
    if t1 - t0 < MIN_COMMON_S:
        return None
    out = Stream([t.copy().trim(t0, t1, nearest_sample=True) for t in pick])
    n = min(t.stats.npts for t in out)
    for t in out:
        t.data = t.data[:n]
    return out


def _write_sac(tr, path, ev, coords, kind):
    from obspy.io.sac import SACTrace
    tr.data = tr.data.astype(np.float32)
    sac = SACTrace.from_obspy_trace(tr)
    sac.reftime = ev["origin"]                       # keeps sample times; b becomes origin-relative
    sac.o = 0.0
    sac.stla, sac.stlo, sac.stel = coords
    sac.evla, sac.evlo, sac.evdp = ev["lat"], ev["lon"], ev["depth"]
    sac.mag = ev["mag"] if np.isfinite(ev["mag"]) else None
    c = tr.stats.channel[2]
    sac.cmpaz, sac.cmpinc = {"Z": (0.0, 0.0), "N": (0.0, 90.0), "E": (90.0, 90.0)}[c]
    sac.idep = "ivel" if kind == "v" else "iacc"
    sac.lcalda = True
    sac.write(path)


class FDSNFetcher:
    """Callable for PocketQuake's orchestrator: fetcher(catalog_csv, waveforms_dir)."""
    backend_name = "fdsn"

    def __init__(self, service: str = "NCEDC", pre_s: float = 240.0, post_s: float = 360.0,
                 workers: int = 2, inventory_xml: str | None = None):
        self.service, self.pre_s, self.post_s = service, pre_s, post_s
        self.workers, self.inventory_xml = workers, inventory_xml

    # ------------------------------------------------------------------ one event
    def _fetch_event(self, inv, sensors, ev, waveforms_dir) -> dict:
        from obspy.clients.fdsn import Client
        from obspy.clients.fdsn.header import FDSNNoDataException
        client = Client(self.service, timeout=600, _discover_services=False)   # one per thread
        eid, origin = ev["event_id"], ev["origin"]
        t0, t1 = origin - self.pre_s, origin + self.post_s
        st = Stream()
        for i in range(0, len(sensors), BULK_CHUNK):
            lines = [(n, s, l or "--", f"{b}?", t0, t1) for n, s, l, b in sensors[i:i + BULK_CHUNK]]
            for attempt in range(3):
                try:
                    st += client.get_waveforms_bulk(lines)
                    break
                except FDSNNoDataException:
                    break
                except Exception:  # noqa: BLE001  (timeouts / 5xx: retry, then fail the event)
                    if attempt == 2:
                        raise
                    _time.sleep(15 * (attempt + 1))
        stage = os.path.join(waveforms_dir, ".staging", eid)
        shutil.rmtree(stage, ignore_errors=True)
        written, missing = [], []
        for n, s, l, b in sensors:
            st3 = st.select(network=n, station=s, location=l, channel=f"{b}?")
            if not st3:
                missing.append(f"{n}.{s}.{l}.{b}")
                continue
            st3.merge(method=1, fill_value=None)
            st3 = _common_stretch(st3, origin)
            if st3 is None:
                missing.append(f"{n}.{s}.{l}.{b}")
                continue
            try:
                for tr in st3:
                    tr.data = tr.data.astype(np.float64)
                    tr.remove_sensitivity(inventory=inv)
                # whatever three components the sensor has (ZNE, Z12, 123, Z23, ...): rotation
                # uses each channel's own azimuth and dip
                st3.rotate("->ZNE", inventory=inv,
                           components=["".join(tr.stats.channel[2] for tr in st3)])
            except Exception as e:  # noqa: BLE001  (metadata gap for this epoch)
                missing.append(f"{n}.{s}.{l}.{b} ({e})")
                continue
            if sorted(tr.stats.channel[2] for tr in st3) != ["E", "N", "Z"]:
                missing.append(f"{n}.{s}.{l}.{b} (not rotatable to ZNE: "
                               f"{','.join(tr.stats.channel for tr in st3)})")
                continue
            kind = "a" if b[1:2] in ("N", "G") else "v"
            z = inv.select(network=n, station=s, location=l, channel=f"{b}?",
                           time=origin)[0][0][0]
            coords = (z.latitude, z.longitude, z.elevation)
            d = os.path.join(stage, f"FDSN.{kind}", "SAC", b)
            os.makedirs(d, exist_ok=True)
            for tr in st3:
                ts = tr.stats.starttime
                name = (f"{n}.{s}.{l}.{tr.stats.channel}.D.{ts.year}."
                        f"{ts.julday:03d}.{ts.strftime('%H%M%S')}.SAC")
                _write_sac(tr, os.path.join(d, name), ev, coords, kind)
            written.append(f"{n}.{s}.{l}.{b}")
        final = os.path.join(waveforms_dir, eid)
        os.makedirs(final, exist_ok=True)
        for kind in ("v", "a"):
            dst = os.path.join(final, f"FDSN.{kind}")
            shutil.rmtree(dst, ignore_errors=True)                 # an interrupted earlier try
            src = os.path.join(stage, f"FDSN.{kind}")
            if os.path.isdir(src):
                os.replace(src, dst)
        shutil.rmtree(stage, ignore_errors=True)
        rec = dict(event_id=eid, origin=str(origin), window_s=[-self.pre_s, self.post_s],
                   service=self.service, sensors=len(sensors), written=len(written),
                   missing=missing)
        with open(os.path.join(final, DONE), "w") as fh:
            json.dump(rec, fh, indent=1)
        return rec

    # --------------------------------------------------------------- all events
    def __call__(self, catalog_csv: str, waveforms_dir: str) -> dict:
        inv_path = self.inventory_xml or os.path.join(os.path.dirname(os.path.abspath(
            waveforms_dir)), "station_table", "inventory.xml")
        inv = read_inventory(inv_path)
        sensors = _sensors(inv)
        cat = pd.read_csv(catalog_csv)
        evs = [dict(event_id=event_id(r.time), origin=UTCDateTime(r.time), lat=float(r.latitude),
                    lon=float(r.longitude), depth=float(r.depth), mag=float(r.magnitude))
               for r in cat.itertuples()]
        todo = [e for e in evs
                if not os.path.exists(os.path.join(waveforms_dir, e["event_id"], DONE))]
        print(f"[fdsn] {len(evs)} events, {len(evs) - len(todo)} already fetched, {len(todo)} to "
              f"fetch; {len(sensors)} sensors; window -{self.pre_s:g}/+{self.post_s:g} s; "
              f"{self.workers} at a time from {self.service}", flush=True)
        t_start, done, failed = _time.time(), 0, []
        with ThreadPoolExecutor(self.workers) as ex:
            futs = {ex.submit(self._fetch_event, inv, sensors, e, waveforms_dir): e
                    for e in todo}
            for f in as_completed(futs):
                e = futs[f]
                try:
                    rec = f.result()
                    done += 1
                    rate = (_time.time() - t_start) / done
                    eta = UTCDateTime() + rate * (len(todo) - done) + 9 * 3600
                    print(f"[fdsn] {e['event_id']} {rec['written']}/{rec['sensors']} sensors | "
                          f"{done}/{len(todo)} done, {rate:.0f} s/event, ETA "
                          f"{eta.strftime('%m-%d %H:%M')} KST", flush=True)
                except Exception as ex_:  # noqa: BLE001
                    failed.append(e["event_id"])
                    print(f"[fdsn] {e['event_id']} FAILED: {ex_}", flush=True)
        if failed:
            print(f"[fdsn] {len(failed)} event(s) failed (rerun to retry): {failed[:10]}")
        return dict(events=len(evs), fetched=done, failed=failed)
