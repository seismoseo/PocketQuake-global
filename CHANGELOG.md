# Changelog

Versioning follows [Semantic Versioning](https://semver.org/). The single source of truth is
`pqglobal.__version__`; `pyproject.toml` reads it.

## 0.3.0 — 2026-11-28

**Benchmark: the November 2026 Pinnacles sequence** (`examples/pinnacles2026.sh`): the 176 USGS ComCat events
within 15 km of the Mw 4.98 mainshock of 20 Nov 2026 on the creeping section of the San Andreas fault, end
to end, with `examples/pinnacles2026/check.py` (the benchmark numbers beside the reference values) and
`examples/pinnacles2026/figure.py` (the README figure). Inputs: the ComCat export as of 2026-11-27, the 10
east-side stations newer than NCSN's alternate list, the USGS NSHM23 fault trace; NCSN's HYPOINVERSE
package is downloaded from NCEDC and checked.

**Added**
- `--velmodel ncsn:<code>` with `--ncsn-multmod` (NCSN's `multmod5a.hyp`) and `--alternate-extra`: NCSN's
  own setup for a region with two models across a fault (e.g. `ncsn:gab`, Bear Valley: Gabilan for the
  Pacific side, Diablo for the stations in the alternate list). HYPOINVERSE reads NCSN's gradient models,
  station delays, P-to-S ratios and alternate list exactly as NCSN does; HypoDD gets the slowness average of
  the two models in at most 20 layers; SKHASH traces takeoff angles through the average and both models.
- `--onset-mag` (default 2.0): every event pair that includes an event of at least this magnitude is
  cross-correlated on the onset only (+-0.05 s around the pick, 1-40 Hz, instead of +-0.5 s, 5-20 Hz), so
  the larger events are placed by where their rupture started; `none` turns it off.
- Offline tests for both.

**Changed**
- The onset window is on by default. Runs made with 0.2.0 or earlier (such as the San Ramon validation
  below) did not use it; pass `--onset-mag none` to repeat them exactly.

**Depends on** PocketQuake v1.18.0 (pipeline: several HYPOINVERSE models per run, onset window for global
clusters, `XCORR_VRAM_FRACTION` / `XCORR_VRAM_CAP` for a shared GPU).

## 0.2.0 — 2026-11-22

- Forward `--s-on-t`, `--orientation-csv`, `--orientation-missing {error,zero,skip}`, and
  `--rotation-reference LAT,LON` to PocketQuake's orchestrator. A table enables S
  cross-correlation on T; the default remains S on the better of N/E, with P on Z.
- Document that orientation tables for FDSN inputs correct residual misorientation after
  metadata-based Z/N/E rotation, to avoid applying metadata orientation twice.
- Add offline CLI forwarding and reference-coordinate validation tests.

**Depends on** PocketQuake v1.17.0 for transverse-component processing.

## 0.1.0 — 2026-10-23

First release: PocketQuake's chain for any region with FDSN data.

**Added**
- `pqglobal.catalog`: USGS ComCat queries (box or circle, CSV format) and ComCat CSV files, written
  as the pipeline's UTC catalog; of two events in the same UTC second the smaller is dropped and
  listed.
- `pqglobal.network_picks`: the network's own picks for the catalog's events (preferred-origin
  arrivals from the FDSN event service, one request per day, cached per day).
- `pqglobal.stations`: the station rule (stations the network picked, within a radius, with a
  three-component sensor at >= 100 samples/s; one sensor per station), with `sensors.csv` giving
  each picked sensor's verdict and `inventory.xml` for the fetcher.
- `pqglobal.waveforms.FDSNFetcher`: event waveforms by FDSN dataselect (origin -240 s .. +360 s),
  signed-sensitivity correction, rotation to Z/N/E from the metadata, SAC in PocketQuake's event
  layout; finished events are skipped on a rerun.
- `pqglobal.velmodel`: gradient models cut into HYPOINVERSE's at most 20 layers, with the
  travel-time error against the gradient model; NCSN CON built in.
- `pqglobal.run` / `pocketquake-global`: catalog → picks → stations → PocketQuake's orchestrator,
  dry run by default; `--velmodel` is required (`con` or a CSV of layers).
- `examples/sanramon2015.sh`, offline tests.

**Depends on** PocketQuake v1.16.0 (submodule `external/PocketQuake`), whose pipeline is
korea-cluster-relocation d067985.

**Verified**: the 30 largest events of the 2015 San Ramon swarm, end to end (picks, HYPOINVERSE,
HypoDD dt.ct and dt.cc, SKHASH, bootstrap errors, notebook, PDF report). Against the NCSS analyst
picks: P 96 % / S 97 % found within 0.5 s at our stations. Against LDEO's NCAeqDD catalog: cluster
shape within 147 m horizontally and 123 m vertically (median), level with SeisComP's rtDD on the
same events (133 / 125 m).
