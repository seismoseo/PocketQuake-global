# Changelog

Versioning follows [Semantic Versioning](https://semver.org/). The single source of truth is
`pqglobal.__version__`; `pyproject.toml` reads it.

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
