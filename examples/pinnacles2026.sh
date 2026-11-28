#!/bin/bash
# The November 2026 Pinnacles sequence (central California), the main benchmark: the 176 USGS ComCat events
# within 15 km of the Mw 4.98 mainshock of 20 Nov 2026 (nc75468702), 1-28 Nov 2026, on the creeping section of
# the San Andreas fault. NCEDC picks, stations and waveforms; NCSN's own Bear Valley velocity setup.
#
#   ./examples/pinnacles2026.sh            dry run: catalog, network picks, station set, then the plan
#   ./examples/pinnacles2026.sh apply      + download waveforms and run the whole pipeline
#
# Inputs (examples/pinnacles2026/README.md): the ComCat catalog as exported on 2026-11-27 (kept, so the
# benchmark does not change when ComCat does) and the 10 stations east of the fault that are newer than NCSN's
# alternate-model list. NCSN's HYPOINVERSE package (velocity models, station delays, alternate list) is
# downloaded from NCEDC once and checked against its SHA-256.
# On a GPU shared with other jobs, limit its memory, e.g.  XCORR_VRAM_FRACTION=0.1 XCORR_VRAM_CAP=0.15
set -euo pipefail

NAME=pinnacles2026              # cluster name
REGION=Pinnacles2026            # display name / source folder prefix
DATA=examples/pinnacles2026
CATALOG=$DATA/comcat_sequence.csv
RADIUS_KM=60                    # stations the network picked, within this distance of the catalog centroid
MIN_SPS=100                     # three-component sensors at this sampling rate or more
VELMODEL=ncsn:gab               # NCSN Bear Valley: Gabilan (Pacific side) + Diablo for the alternate list
ALT_EXTRA=$DATA/alternate_extra_newer_stations.txt
ONSET_MAG=2.0                   # pairs with an event of M >= 2 are cross-correlated on the onset only
CORES=${CORES:-2}               # cross-correlation worker cap (memory on a shared machine)
NCSN_URL=https://ncedc.org/pub/catalogs/ncsn/hypoinverse/hypoinv-standalone-20230210.zip
NCSN_SHA256=863215028a53bd9dc2a285498890e8a6cdee187f8937d9e8d94def63a3271234

cd "$(dirname "$0")/.."
WORK=work/$NAME
NCSN=work/ncsn_hypoinverse
if [ ! -f "$NCSN/multmod5a.hyp" ]; then
    mkdir -p "$NCSN"
    curl -fsSL -o "$NCSN/hypoinv-standalone.zip" "$NCSN_URL"
    if command -v sha256sum > /dev/null; then sum=$(sha256sum "$NCSN/hypoinv-standalone.zip")
    else sum=$(shasum -a 256 "$NCSN/hypoinv-standalone.zip"); fi
    [ "${sum%% *}" = "$NCSN_SHA256" ] || { echo "NCSN package checksum mismatch: $NCSN_URL"; exit 1; }
    unzip -q -o "$NCSN/hypoinv-standalone.zip" -d "$NCSN"
fi
ARGS=(--name "$NAME" --region "$REGION" --work "$WORK" --catalog-csv "$CATALOG"
      --radius-km "$RADIUS_KM" --min-sps "$MIN_SPS" --cores "$CORES"
      --velmodel "$VELMODEL" --ncsn-multmod "$NCSN/multmod5a.hyp" --alternate-extra "$ALT_EXTRA"
      --onset-mag "$ONSET_MAG")

if [ "${1:-}" = "apply" ]; then
    python -m pqglobal.run "${ARGS[@]}" --apply
else
    python -m pqglobal.run "${ARGS[@]}" --stop-after stations
    python -m pqglobal.run "${ARGS[@]}"           # prints the cluster settings; nothing is downloaded
    echo "dry run done: inspect $WORK/station_table/sensors.csv, then rerun with 'apply'"
fi
