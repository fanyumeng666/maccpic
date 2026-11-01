#!/bin/sh
# T8: a straight run of 6 periods and a run stopped after 3 periods + restart to 6 must
#     give bit-identical history (except the wall-clock column), profiles,
#     distributions and x-t maps.  2 ranks, debug executable, n_avg = 4 so
#     that the measurement window is open when the restart file is written.
cd "$(dirname "$0")" || exit 1
BIN=$(cd ../.. && pwd)/ccp1d_dbg
MPIRUN="mpirun --allow-run-as-root --oversubscribe -np 2"
rm -rf straight split; mkdir straight split
cat > straight/ccp1d.in <<INP
&ccp1d
  L = 0.067d0, nx = 32, freq = 13.56d6, V0 = 450d0, n_gas = 0d0
  n0 = 2.56d14, Te0 = 30000d0, Ti0 = 300d0, nppc = 64
  steps_per_period = 200, n_periods = 6, n_avg = 4, n_phase = 50, seed = 4321, restart_every = 3
/
INP
sed 's/n_periods = 6/n_periods = 6, periods_this_run = 3/' straight/ccp1d.in > split/ccp1d.in
(cd straight && $MPIRUN $BIN > run.log 2>&1) || { echo "    FAIL T8 straight run failed"; exit 1; }
(cd split && $MPIRUN $BIN > run1.log 2>&1) || { echo "    FAIL T8 first segment failed"; exit 1; }
sed -i 's/periods_this_run = 3/periods_this_run = 0, restart = .true./' split/ccp1d.in
(cd split && $MPIRUN $BIN > run2.log 2>&1) || { echo "    FAIL T8 restart segment failed"; exit 1; }
awk '!/^#/{$26=""; print}' straight/history.dat > straight/h.txt
awk '!/^#/{$26=""; print}' split/history.dat > split/h.txt
fail=0
if cmp -s straight/h.txt split/h.txt; then echo "    PASS T8 history identical (7 lines)"; else echo "    FAIL T8 history differs"; diff straight/h.txt split/h.txt | head -6; fail=1; fi
for f in profiles.dat iedf_left.dat iedf_right.dat eedf_left.dat eepf.dat xt_phi.dat xt_ne.dat xt_Je.dat xt_PeE.dat; do
  if cmp -s straight/$f split/$f; then echo "    PASS T8 $f identical"; else echo "    FAIL T8 $f differs"; fail=1; fi
done
if [ $fail -eq 0 ]; then echo "test_restart: PASS"; else echo "test_restart: FAIL"; exit 1; fi
