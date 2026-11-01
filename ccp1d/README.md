# ccp1d — 1D bounded electrostatic PIC/MCC for an RF capacitive discharge

Fortran 90 + MPI. Target: the helium benchmark of Turner et al., Phys. Plasmas 20, 013507 (2013), Case 1.
Specification: [`../HANDOVER.md`](../HANDOVER.md); plan and test matrix: [`PLAN.md`](PLAN.md); progress: [`PROGRESS.md`](PROGRESS.md).

**Status (2026-11-01): M0 + M1 done** — collisionless code (push, CIC deposit, Poisson, absorbing
electrodes, diagnostics, restart) with passing unit tests and the collisionless-sheath test.
Collisions (M2, `xsec.f90` + `mcc.f90`) are next; running Case 1 now gives a collisionless discharge.

## Build, test, run

```sh
./setup.sh          # once per fresh container (repository root): gfortran, Open MPI, numpy, matplotlib
cd ccp1d
make                # -O3 executable ./ccp1d (objects in obj/)
make debug          # checked executable ./ccp1d_dbg (-fcheck=all -ffpe-trap ...)
make check          # unit tests T1-T3, T8, T10 with the debug flags (~1 min)
make check-rng      # statistical tests of the random generator in the loading pattern (~1 min)
mkdir -p runs/case1_test && cp inputs/case1_test.in runs/case1_test/ccp1d.in
cd runs/case1_test && mpirun --allow-run-as-root -np 2 ../../ccp1d [input_file]
```

The run directory holds `ccp1d.in` (namelist `&ccp1d`, see `inputs/*.in` and `input.f90` for every key)
and all output. Restart: `restart = .true.` continues from `restart_*.bin`; `periods_this_run = N`
stops after N periods, `max_wall_minutes = M` stops cleanly after M minutes (both write restart files).

## Source files

| file | role |
| --- | --- |
| `input.f90` | namelist, checks, derived quantities (dx, dt, weight, maxnp), drive voltage `v_drive(t)` |
| `global.f90` | constants, species, grid and particle arrays, counters, `add/remove_particle`, `abort_run` |
| `random.f90` | ran2 (as em2's dran2) with saveable state, Box-Muller |
| `poisson.f90` | tridiagonal Poisson solve with Dirichlet electrodes, E on nodes and faces |
| `deposit.f90` | CIC deposit of rho (all ranks) and of the species moments |
| `push.f90` | leapfrog push, electrode absorption, exact work and Shockley-Ramo current sums |
| `diag.f90` | history, window averages, histograms, x-t maps, energy balance, info.txt, restart files |
| `init.f90` | random uniform loading with Maxwellian velocities, field at t = 0, leapfrog half step |
| `main.f90` | MPI, time loop |
| `tests/` | T1 poisson, T2 deposit, T3 random, T8 restart (bit-identical), T10 Langmuir oscillation; `random/test_tails`, `test_seeds`, `test_corr` (generator statistics, `make check-rng`) |
| `tools/` | `plot_history.py`, `check_sheath.py`, `plot_initial.py`, `common.py` |
| `inputs/` | `case1.in`, `case1_test.in`, `sheath.in`, `init_case1.in` |
| `doc/` | `initial_conditions.md` (what is loaded at t = 0 and why, with figures) |
| `results/` | small result files of the reference runs (`sheath/`, `init_case1/`) |

## Output files (all ASCII, `#` headers with column names and units)

`info.txt` (parameters, derived numbers, results and the energy balance), `history.dat` (one line per period),
`profiles.dat` (time averages over the last `n_avg` periods), `iedf_left/right.dat`, `eedf_left/right.dat`,
`eepf.dat`, `xt_*.dat` (phase-resolved maps), `restart_*.bin`. Column lists are in the file headers and in
`tools/common.py`.
