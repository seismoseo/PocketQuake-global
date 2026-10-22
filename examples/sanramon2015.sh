#!/bin/bash
# San Ramon (California) 2015 swarm, the validation case: 680 USGS ComCat events in the box of the
# SeisComP/rtDD study, Oct 1 - Dec 1 2015, NCEDC stations and waveforms, NCSN CON velocity model.
#
#   ./examples/sanramon2015.sh            dry run: catalog, network picks, station set, then the plan
#   ./examples/sanramon2015.sh apply      + download waveforms and run the whole pipeline
#   LARGEST=30 ./examples/sanramon2015.sh apply   a quick subset: the 30 largest events
set -euo pipefail

NAME=sanramon2015               # cluster name (subsets get a suffix below)
REGION=SanRamon2015             # display name / source folder prefix
BOX=37.75,37.83,-122.01,-121.91 # lat0,lat1,lon0,lon1 of the ComCat query
START=2015-10-01
END=2015-12-01
RADIUS_KM=60                    # stations the network picked, within this distance of the catalog centroid
MIN_SPS=100                     # three-component sensors at this sampling rate or more
VELMODEL=con                    # NCSN Concord-Calaveras (built in) or a CSV of layers: top_km,vp_km_s
VPVS=1.747
CORES=2                        # cross-correlation worker cap (memory on a shared machine)
LARGEST=${LARGEST:-}            # e.g. 30 for a quick subset

cd "$(dirname "$0")/.."
WORK=work/$NAME${LARGEST:+_top$LARGEST}
ARGS=(--name "$NAME${LARGEST:+top$LARGEST}" --region "$REGION${LARGEST:+Top$LARGEST}" --work "$WORK"
      --usgs-box "$BOX" --start "$START" --end "$END" --radius-km "$RADIUS_KM" --min-sps "$MIN_SPS"
      --velmodel "$VELMODEL" --vpvs "$VPVS" --cores "$CORES")
[ -n "$LARGEST" ] && ARGS+=(--largest "$LARGEST")

if [ "${1:-}" = "apply" ]; then
    python -m pqglobal.run "${ARGS[@]}" --apply
else
    python -m pqglobal.run "${ARGS[@]}" --stop-after stations
    python -m pqglobal.run "${ARGS[@]}"           # prints the cluster settings; nothing is downloaded
    echo "dry run done: inspect $WORK/station_table/sensors.csv, then rerun with 'apply'"
fi
