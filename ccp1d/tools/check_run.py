#!/usr/bin/env python3
"""check_run.py <run_dir> [--full]: sanity checks of a benchmark run (T13; the case is identified from n_gas and V0).
Default tolerances are for the test mode (nx = 32); --full uses the full-size ones."""
import os
import sys
import numpy as np
from common import load_dat, load_info, history_columns, profile_columns, check, check_range, expected_phase_drive, turner_table3, QE

args = [a for a in sys.argv[1:] if not a.startswith('--')]
full = '--full' in sys.argv
d = args[0] if args else '.'
h, _ = load_dat(os.path.join(d, 'history.dat'))
c = history_columns()
prof, _ = load_dat(os.path.join(d, 'profiles.dat'))
pc = profile_columns()
info = load_info(os.path.join(d, 'info.txt'))
ok = True
tol_drift, tol_nc, tol_asym, tol_flux, tol_clos, tol_pow = (0.02, 0.05, 0.02, 0.05, 1e-3, 0.05) if full else (0.20, 1.0, 0.05, 0.10, 5e-3, 0.10)

# 1 finite
bad = 0
for f in ['history.dat', 'profiles.dat', 'iedf_left.dat', 'eepf.dat']:
    a, _ = load_dat(os.path.join(d, f))
    bad += int(np.sum(~np.isfinite(a)))
ok &= check('1  non-finite entries in the output files', bad, 0)
# 2 steady state: drift of N_e, N_i, n_i_centre over the last 20 % of the periods (per 100 periods)
n = len(h); k = max(3, n // 5)
p = h[-k:, c['period']]
worst = 0
for name in ['N_e', 'N_i', 'n_i_center']:
    y = h[-k:, c[name]]
    slope = np.polyfit(p, y, 1)[0]
    worst = max(worst, abs(slope * 100 / max(y.mean(), 1e-300)))
ok &= check('2  max drift of N_e, N_i, n_i_c over the last 20 % (fraction per 100 periods)', worst, tol_drift)
# 3 charge imbalance: sheaths make N_i > N_e, by a few to ~15 %
r = (h[-1, c['N_i']] - h[-1, c['N_e']]) / max(h[-1, c['N_i']], 1.0)
if h[-1, c['N_e']] == 0:
    print('    FAIL 0  the plasma is extinguished (N_e = 0 at the end): the remaining checks are not meaningful')
    ok = False
ok &= check_range('3  (N_i - N_e)/N_i (positive: the sheaths hold ions only)', r, 0.0, 0.6)
# 4 phi(0, t) = V0 sin(2 pi f t) from xt_phi
xt, hdr = load_dat(os.path.join(d, 'xt_phi.dat'))
nph = xt.shape[0]
V0 = info['V0']
Vmax = V0 + abs(info.get('V0_lf', 0.0))          # largest electrode voltage (dual frequency: both terms)
expect = expected_phase_drive(info, nph)
err = np.max(np.abs(xt[:, 0] - expect)) / Vmax
ok &= check('4  phi(0,t) vs the drive averaged over the steps of each phase bin, relative', err, 1e-6)   # files carry 7 digits
# 5 centre density vs Turner: Case 1 from the reference profile, Cases 2-4 from Table III
x = prof[:, pc['x']]; ni = prof[:, pc['n_i']]; ne = prof[:, pc['n_e']]
mid = len(x) // 2
nc = ni[mid - 2:mid + 3].mean()
case, row = turner_table3(info)
if case == 0:
    print('       5  no benchmark case matches n_gas = %.3e, V0 = %.0f V: centre n_i = %.3e (information only)' % (info['n_gas'], info['V0'], nc))
else:
    nref = 1.405e14 if case == 1 else row[3]
    ok &= check('5  |n_i centre / Turner case %d (%.3e) - 1|' % (case, nref), abs(nc / nref - 1), tol_nc)
# 6 symmetry (half-domain ion numbers)
nl, nr = np.trapezoid(ni[:mid + 1], x[:mid + 1]), np.trapezoid(ni[mid:], x[mid:])
ok &= check('6  |N_i,left - N_i,right|/(N_i,left + N_i,right)', abs(nl - nr) / (nl + nr), tol_asym)
# 7 flux balance per electrode and left/right
for side in ['left', 'right']:
    ok &= check('7  |Gamma_e/Gamma_i - 1| ' + side, abs(info['Gamma_e_' + side] / max(info['Gamma_i_' + side], 1e-300) - 1), tol_flux)
ok &= check('7  |Gamma_i,left/Gamma_i,right - 1|', abs(info['Gamma_i_left'] / max(info['Gamma_i_right'], 1e-300) - 1), tol_flux)
# 8 ion impact energies
iedf, _ = load_dat(os.path.join(d, 'iedf_left.dat'))
e_bin, f = iedf[:, 0], iedf[:, 1]
emax = e_bin[np.nonzero(f)[0][-1]] if np.any(f > 0) else 0.0
ok &= check('8  max ion impact energy - (V0 + V0_lf) [eV]', emax - Vmax, 50.0)
ok &= check_range('8  mean ion impact energy / (V0 + V0_lf) (left)', info.get('E_i_mean_left', 0.0) / Vmax, 0.1, 0.6)
# 9 electron energy at the centre
ok &= check_range('9  <eps_e> centre [eV]', info['eps_e_center'], 3.0, 20.0)
# 10 energy bookkeeping and the three powers
ok &= check('10 closure_relative (bookkeeping identity)', info['closure_relative'], tol_clos)
ptot = info['P_part_total'] if info['P_part_total'] != 0 else 1e-300
ok &= check('10 |P_work integral / P_part - 1| (exact work deposit; info.txt has 7 digits)', abs((info['P_work_e_integral'] + info['P_work_i_integral']) / ptot - 1), 2e-6)
if info['P_part_e'] != 0 and info['P_part_i'] != 0:
    print('       grid J.E: P_JE_e/P_part_e = %.3f, P_JE_i/P_part_i = %.3f (centering error of the grid product, information only)' %
          (info['P_JE_e'] / info['P_part_e'], info['P_JE_i'] / info['P_part_i']))
ok &= check('10 |P_VI/P_part - 1|', abs(info['P_VI'] / ptot - 1), tol_pow)
ok &= check('10 P_numerical_relative (energy error of the scheme)', info['P_numerical_relative'], tol_pow)
# 11 null-collision guard
ok &= check('11 numax exceeded count', info['numax_exceeded_count'], 0)
# 12 ionization balance: ionizations per period vs ion losses per period
loss = (info['Gamma_i_left'] + info['Gamma_i_right']) * info['period'] / info['weight']
ok &= check('12 |ionizations per period / ion losses per period - 1|', abs(info['ionizations_per_period'] / max(loss, 1e-300) - 1), tol_flux)
# 13 timing
print('    13 wall time: %.1f s for %d periods (%.3f s/period), %d ranks' %
      (h[-1, c['wall']], int(h[-1, c['period']]), h[-1, c['wall']] / max(1, h[-1, c['period']]), info['nprocs']))
print('check_run (%s):' % ('full-size tolerances' if full else 'test-mode tolerances'), 'PASS' if ok else 'FAIL')
raise SystemExit(0 if ok else 1)
