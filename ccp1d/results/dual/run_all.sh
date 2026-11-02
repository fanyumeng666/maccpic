#!/bin/sh
# M7 dual-frequency runs, one after the other (4 ranks each)
cd "$(dirname "$0")" || exit 1
BIN=$(cd ../.. && pwd)/ccp1d
MPIRUN="mpirun --allow-run-as-root --oversubscribe -np 4"
for r in dual_100 dual_200 dual_300 dual_200_lf; do
  mkdir -p $r && sed 's#../../data/he#../../../data/he#' ../../inputs/$r.in > $r/ccp1d.in
  echo "==== $r start $(date +%H:%M:%S)"
  (cd $r && $MPIRUN $BIN > run.log 2>&1)
  echo "==== $r done $(date +%H:%M:%S): $(tail -1 $r/run.log)"
done
echo "ALL DUAL RUNS DONE"
