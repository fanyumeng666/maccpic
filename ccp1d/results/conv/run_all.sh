#!/bin/sh
# M5 convergence runs, one after the other (4 ranks each): dt/2, dx/2, other seed, nppc 256, nppc 1024
cd "$(dirname "$0")" || exit 1
BIN=$(cd ../.. && pwd)/ccp1d
MPIRUN="mpirun --allow-run-as-root --oversubscribe -np 4"
make_in() {  # name, sed expression applied to inputs/case1.in
  mkdir -p "$1"; sed -e "$2" -e 's/max_wall_minutes = 100d0/max_wall_minutes = 100d0, xsec_dir = "..\/..\/..\/data\/he"/' ../../inputs/case1.in > "$1/ccp1d.in"
}
make_in dt2     's/steps_per_period = 400/steps_per_period = 800/'
make_in nx256   's/nx = 128/nx = 256/'
make_in seed2   's/seed = 12345/seed = 777/'
make_in nppc256 's/nppc = 512/nppc = 256/'
make_in nppc1024 's/nppc = 512/nppc = 1024/'
for r in dt2 nx256 seed2 nppc256 nppc1024; do
  echo "==== $r start $(date +%H:%M:%S)"
  (cd $r && $MPIRUN $BIN > run.log 2>&1)
  echo "==== $r done $(date +%H:%M:%S): $(tail -1 $r/run.log)"
done
echo "ALL CONVERGENCE RUNS DONE"
