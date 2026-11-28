# The Pinnacles 2026 benchmark

The main benchmark of PocketQuake-global: the 176 USGS ComCat events within 15 km of the Mw 4.98
mainshock of 20 November 2026 (nc75468702), on the creeping section of the San Andreas fault near
Pinnacles, California, relocated end to end and checked against NCSN's analyst picks, NCSN's catalog, the
USGS fault trace and LDEO's real-time double-difference catalog. The results are in the main
[README](../../README.md#benchmark-the-2026-pinnacles-sequence); what you need before you run it (GPU,
about 170 GB of disk, network access to NCEDC) is listed there too.

## Run it

From the repository root, with the environment and tools of the main README's Install section:

```bash
./examples/pinnacles2026.sh                    # dry run: catalog, NCSN picks, the 42 stations, the plan
CORES=6 ./examples/pinnacles2026.sh apply      # download and the whole chain
python examples/pinnacles2026/check.py         # this run's numbers beside the reference values
python examples/pinnacles2026/figure.py        # work/pinnacles2026/pinnacles2026.png (PyGMT)
python examples/pinnacles2026/compare_ddrt.py  # work/pinnacles2026/ldeo_ddrt/ (numbers and figures, PyGMT)
```

`ALT_EXTRA= ./examples/pinnacles2026.sh apply` runs the variant with NCSN's alternate list as distributed
(cluster `pinnacles2026_ncsnlist`); point `check.py` at it with `--run-dir` and `--work`.

## The files

| file | what it is |
|---|---|
| [`../pinnacles2026.sh`](../pinnacles2026.sh) | the run: every setting of the benchmark in one place |
| `comcat_sequence.csv` (+ `.request.json`) | the catalog: ComCat's export, kept unchanged |
| `alternate_extra_newer_stations.txt` | the 10 east-side stations newer than NCSN's alternate list |
| `usgs_nshm23_san_andreas.geojson` (+ `.request.json`) | the San Andreas fault trace, USGS 2023 National Seismic Hazard Model fault sections |
| `ldeo_ddrt_sequence.txt` (+ `.request.json`) | the reference relocation: LDEO's real-time double-difference catalog for these events |
| `reference.json` | the reference values that `check.py` prints beside a run's own |
| `check.py` | the benchmark numbers of a run: picks against NCSN's analyst picks, located events, mechanisms, the line of epicenters, the comparison with LDEO |
| `figure.py` | the README's header figure: stations by velocity model, NCSN's catalog against the relocation |
| `compare_ddrt.py` | NCSN's catalog, LDEO's catalog and the relocation on the same events: offset, relative differences, width of the line, the fault plane of the aftershocks; two figures |

Each `.request.json` records where its file came from (URL or query), when it was retrieved, and its
SHA-256, so a copy can be checked against the original.

## The inputs

**`comcat_sequence.csv`**: the USGS ComCat events from 1 to 27 November 2026 (query 1 November 00:00 to
28 November 00:00 UTC) within 15 km of the Mw 4.98 mainshock (nc75468702, 2026-11-20 06:49:16.34 UTC,
36.6072° N, 121.2083° W): 176 events, M 0.7–4.98. This is the CSV exactly as ComCat returned it on
2026-11-27 14:50 UTC. The file is kept so the benchmark stays the same when ComCat revises the catalog.

**`alternate_extra_newer_stations.txt`**: 10 stations, in NCSN's alternate-list format (`STA   NET  A`).
In Bear Valley, NCSN locates events with the Gabilan model (GAB) for stations on the Pacific side of the
San Andreas fault, and with the Diablo model (DIA) for the stations in its alternate list
(`hypoinverse/models/multmod5.alt`), the stations east of the fault. NCSN's file says it was updated "for
stations installed up to spring 2013", so it lacks the east-side stations installed since. Without them,
these stations get the Pacific-side model, and their P residuals arrive up to 1.3 s late.

These 10 are the stations of this benchmark that lie east of the fault and are missing from NCSN's list.
The side is that of the nearest segment of the USGS NSHM23 San Andreas fault trace.

| station | distance east of the fault (km) | installed |
|---|---:|---|
| NC.BPV | 4.0 | 2019-03-13 |
| NC.HBH | 5.8 | 2019-03-13 |
| NC.HHP | 6.3 | 2023-04-17 |
| NC.BPA | 10.7 | 2024-10-30 |
| BK.BABI | 17.2 | 2019-05-23 |
| NC.HLO | 18.3 | 2019-05-21 |
| BK.STRA | 19.6 | 2025-11-05 |
| BK.AZUL | 29.9 | 2024-04-02 |
| BK.TUMY | 31.0 | 2024-11-06 |
| WR.CK14 | 47.5 | 2022-01-26 |

With them on DIA, the median HYPOINVERSE RMS falls from 0.28 to 0.19 s, and the relocated sequence lies
0.58 km closer to the mapped fault trace (0.95 km south-west of it), with the same shape.

**`usgs_nshm23_san_andreas.geojson`**: the San Andreas fault sections of the USGS 2023 National Seismic
Hazard Model between 35.4° and 37.6° N, as returned by the USGS query in its request file. `check.py`,
`figure.py` and `compare_ddrt.py` measure distances from this trace.

**`ldeo_ddrt_sequence.txt`**: the reference relocation, from LDEO's real-time double-difference catalog
(Waldhauser 2009). LDEO relocates every NCSN event in real time, from NCSN's picks and waveform
cross-correlation, relative to the events of its 1984–2021 base catalog; its picks, stations, velocity model
and reference events all differ from PocketQuake-global's. The file keeps LDEO's header verbatim, then the
lines of the 176 benchmark events (173 with a double-difference solution, 3 at NCSN's location), from the
catalog file created 2026-11-28 05:50:47 EST. LDEO revises its catalog as NCSN revises events: on that date,
149 of the 152 events after the mainshock were still relocated from NCSN's first, preliminary solution. The
copy is kept so the reference stays the same. Please cite Waldhauser (2009) and the catalog version when
you use it.

**NCSN's HYPOINVERSE package** (`hypoinv-standalone-20230210.zip` from NCEDC; the script downloads it into
`work/ncsn_hypoinverse/` and checks its SHA-256): the gradient models `gab.crt` and `dia.crt`, the station
delays `gab5.del`, the P-to-S ratios, and the alternate list, all read from `multmod5a.hyp` as NCSN uses
them.

## Reading `check.py`

One line per quantity: its value in this run and in the reference run (`reference.json`, one run from an
empty folder on 2026-11-28 with PyTorch 2.11 and CUDA 12.8 on an NVIDIA RTX PRO 6000 Blackwell).

| group | quantities |
|---|---|
| `picks_P`, `picks_S` | PhaseNet+ against NCSN's analyst picks at the same station and phase: pairs, median absolute difference (s), share within 0.1 s |
| `located` | events in the catalog, located by HYPOINVERSE, relocated by HypoDD from catalog (dt.ct) and cross-correlation (dt.cc) differences; HYPOINVERSE's median RMS |
| `mechanisms` | SKHASH: events with a mechanism, of quality A or B (and of those, with a single solution); the mainshock's mechanism and its Kagan angle to the NC moment tensor |
| `shape_dtcc`, `shape_ncsn_catalog` | the line through the epicenters within 6 km of the mainshock: strike, scatter across it (m), the mainshock's distance from it (m; positive south-west, negative north-east), share of aftershocks south-east of the mainshock, median distance south-west of the fault trace (km) |
| `ldeo_ddrt` | against LDEO: events compared, the offset between the two clouds, the median differences left after removing it, then pairs of values as [PocketQuake-global, LDEO]: width of the line (m) and events more than 0.5 km off it; the plane fitted to the aftershocks and its angles to the USGS Mww and NC moment-tensor planes |

Another PyTorch build or GPU can change a few picks, and the numbers with them: counts by a few events,
medians by milliseconds, positions and widths by metres to tens of metres. A rerun in a folder used before
drifts the same way (main README, Known limitations).
