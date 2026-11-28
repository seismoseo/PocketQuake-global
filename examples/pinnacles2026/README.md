# Inputs of the Pinnacles 2026 benchmark

Used by `examples/pinnacles2026.sh`.

**`comcat_sequence.csv`**: USGS ComCat events from 1 to 28 November 2026 within 15 km of the Mw 4.98
mainshock (nc75468702, 2026-11-20 06:49:16.34 UTC, 36.6072° N, 121.2083° W): 176 events, M 0.7–4.98.
This is the CSV exactly as ComCat returned it on 2026-11-27 14:50 UTC. The query, time and SHA-256 are in
`comcat_sequence.request.json`. The file is kept so the benchmark stays the same when ComCat revises
the catalog.

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

With them on DIA, the median HYPOINVERSE RMS falls from 0.28 to 0.19 s, and 56 instead of 48 focal
mechanisms reach quality A or B. The relocated sequence then lies 0.96 km instead of 1.54 km south-west
of the mapped fault trace, with the same shape. Without `--alternate-extra`, the run uses NCSN's list
exactly as distributed.

**NCSN's HYPOINVERSE package** (`hypoinv-standalone-20230210.zip` from NCEDC; the script downloads it and
checks its SHA-256): the gradient models `gab.crt` and `dia.crt`, the station delays `gab5.del`, the P-to-S
ratios, and the alternate list, all read from `multmod5a.hyp` as NCSN uses them.
