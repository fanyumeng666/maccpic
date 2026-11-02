#!/usr/bin/env python3
"""check_diag.py <run_dir>: T16, self-consistency of the diagnostics of a run.
Exact identities between files (same counts written in two places) must hold to
the printed precision; physical consistency checks (total current uniform in x,
sheath edges widest at the right phase, heating at the sheath edge, EEPF slopes
bracketing the mean energy) have loose tolerances.  Writes diag_checks.png."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import load_dat, load_info, profile_columns, run_dir_from_argv, check, check_range, expected_phase_drive, QE, EPS0

d = run_dir_from_argv()
info = load_info(os.path.join(d, 'info.txt'))
prof, _ = load_dat(os.path.join(d, 'profiles.dat'))
pc = profile_columns()
L, V0, T, w, dt = info['L'], info['V0'], info['period'], info['weight'], info['dt']
nx = int(info['nx']); dx = L / nx
ok = True
tol_file = 2e-5      # files carry 7 significant digits; sums of many such numbers


def xt(name):
    a, _ = load_dat(os.path.join(d, 'xt_%s.dat' % name))
    return a


# ---- 1 wall energy distributions integrate to the fluxes in info.txt ----
for f, key in [('iedf_left', 'Gamma_i_left'), ('iedf_right', 'Gamma_i_right'), ('eedf_left', 'Gamma_e_left'), ('eedf_right', 'Gamma_e_right')]:
    a, _ = load_dat(os.path.join(d, f + '.dat'))
    de = a[1, 0] - a[0, 0]
    ok &= check('1  integral of %s / %s - 1' % (f, key), abs(np.sum(a[:, 1]) * de / info[key] - 1), tol_file)
# ---- 2 EEPF normalisation: integral F deps = n_e in the central 10 % ----
ep, _ = load_dat(os.path.join(d, 'eepf.dat'))
de = ep[1, 0] - ep[0, 0]
x = prof[:, pc['x']]
mc = np.abs(x - 0.5 * L) < 0.05 * L
ne_c = np.trapezoid(prof[mc, pc['n_e']], x[mc]) / (x[mc][-1] - x[mc][0])
ok &= check('2  integral F(eps) deps / n_e(central 10 %) - 1', abs(np.sum(ep[:, 1]) * de / ne_c - 1), 0.02)
# EEPF slopes
fp = ep[:, 2]; e = ep[:, 0]
teff = info['T_eff_center']
lo = (e > 1) & (e < 8) & (fp > 0); hi = (e > 25) & (e < 80) & (fp > 0)
t_lo = -1 / np.polyfit(e[lo], np.log(fp[lo]), 1)[0]
t_hi = -1 / np.polyfit(e[hi], np.log(fp[hi]), 1)[0]
print('       EEPF temperatures from the slopes: 1-8 eV: %.2f eV, 25-80 eV: %.2f eV; 2/3 <eps_e> = %.2f eV' % (t_lo, t_hi, teff))
ok &= check_range('2  T_low <= T_eff <= T_high (bi-Maxwellian brackets the mean energy): T_eff/T_low', teff / t_lo, 1.0, 1e9)
ok &= check_range('2  T_eff/T_high', teff / t_hi, 0.0, 1.0)
# ---- 3 x-t maps against profiles.dat and info.txt ----
phi, ne, ni, je, ji, pee, ioniz, pwork = (xt(n) for n in ['phi', 'ne', 'ni', 'Je', 'Ji', 'PeE', 'ioniz', 'Pe_work'])
nph = phi.shape[0]
inner = slice(1, nx)
for name, z, col in [('n_e', ne, 'n_e'), ('n_i', ni, 'n_i'), ('phi', phi, 'phi')]:
    ref = prof[:, pc[col]]
    dev = np.max(np.abs(z.mean(axis=0)[inner] - ref[inner]) / np.max(np.abs(ref)))
    ok &= check('3  phase average of xt_%s = profiles.dat (interior nodes)' % name, dev, tol_file)
# ionizations per period
nion = np.sum(ioniz) * dx * (T / nph) / w
ok &= check('3  ionizations per period from xt_ioniz / info - 1', abs(nion / info['ionizations_per_period'] - 1), tol_file)
# electron power: exact work deposit and grid product
pw = np.mean(np.sum(pwork, axis=1) * dx)
ok &= check('3  <integral P_work_e dx> from xt_Pe_work / P_part_e - 1', abs(pw / info['P_part_e'] - 1), tol_file)
pje = np.mean(dx * (np.sum(pee[:, 1:nx], axis=1) + 0.5 * (pee[:, 0] + pee[:, nx])))
ok &= check('3  <integral J_e.E dx> from xt_PeE / P_JE_e - 1', abs(pje / info['P_JE_e'] - 1), tol_file)
# phase column 0 = V(t) (bin-averaged over the steps)
Vmax = V0 + abs(info.get('V0_lf', 0.0))
expect = expected_phase_drive(info, nph)
ok &= check('3  xt_phi(:,0) = V(t) averaged over the steps of each bin, relative', np.max(np.abs(phi[:, 0] - expect)) / Vmax, 1e-6)
dual = info.get('V0_lf', 0.0) != 0.0
# ---- 4 total current J_e + J_i + eps0 dE/dt uniform in x at every phase ----
E = -np.gradient(phi, dx, axis=1)
dEdt = (np.roll(E, -1, axis=0) - np.roll(E, 1, axis=0)) / (2 * T / nph)
it_ = je + ji + EPS0 * dEdt
amp = np.max(np.abs(it_.mean(axis=1)))
spread = np.sqrt(np.mean((it_[:, inner] - it_[:, inner].mean(axis=1, keepdims=True))**2))
print('       total current: amplitude %.1f A/m^2, rms variation along x %.1f A/m^2' % (amp, spread))
ok &= check('4  rms x-variation of the total current / its amplitude (continuity + Gauss)', spread / amp, 0.15)
# ---- 5 sheath edges: widest when the electrode is most negative ----
t = (np.arange(nph) + 0.5) / nph
r = ne / np.maximum(ni, 1e-300)
sl = np.array([x[np.nonzero(r[p] >= 0.5)[0][0]] for p in range(nph)])
sr = np.array([x[np.nonzero(r[p] >= 0.5)[0][-1]] for p in range(nph)])
print('       sheath edge (n_e/n_i = 0.5): left %.2f-%.2f cm, right %.2f-%.2f cm' % (sl.min() * 1e2, sl.max() * 1e2, sr.min() * 1e2, sr.max() * 1e2))
ok &= check('5  left sheath oscillates (excursion in cells)', -(sl.max() - sl.min()) / dx, -2.0)
ph_l = t[np.argmax(sl)]; ph_r = t[np.argmin(sr)]
if dual:
    print('       dual frequency: widest left sheath at t/T = %.2f, right at %.2f (phase test applies to one frequency only)' % (ph_l, ph_r))
else:
    ok &= check('5  |phase of the widest left sheath - 0.75| (V = -V0 at t/T = 0.75)', min(abs(ph_l - 0.75), 1 - abs(ph_l - 0.75)), 0.15)
    ok &= check('5  |phase of the widest right sheath - 0.25| (plasma potential highest at t/T = 0.25)', min(abs(ph_r - 0.25), 1 - abs(ph_r - 0.25)), 0.15)
# ---- 6 electron heating peaks within the sheath-edge excursion ----
pw_prof = prof[:, pc['P_work_e']]
mid = nx // 2
xpk_l = x[np.argmax(pw_prof[:mid])]; xpk_r = x[mid + np.argmax(pw_prof[mid:])]
print('       phase-averaged electron heating peaks at x = %.2f and %.2f cm' % (xpk_l * 1e2, xpk_r * 1e2))
ok &= check_range('6  left heating peak within the left sheath-edge excursion (0 = min, 1 = max)', (xpk_l - sl.min()) / (sl.max() - sl.min()), -0.2, 1.2)
ok &= check_range('6  right heating peak within the right sheath-edge excursion', (sr.max() - xpk_r) / (sr.max() - sr.min()), -0.2, 1.2)
# ---- 7 ion energies: mean below the time-averaged sheath voltage, max below V0 ----
phic = prof[mid - 2:mid + 3, pc['phi']].mean()
iedf, _ = load_dat(os.path.join(d, 'iedf_left.dat'))
emax = iedf[np.nonzero(iedf[:, 1])[0][-1], 0]
print('       time-averaged plasma potential %.1f V; mean ion impact energy %.1f eV; maximum %.0f eV' % (phic, info['E_i_mean_left'], emax))
ok &= check_range('7  mean ion impact energy / e phi_c (collisional sheath: well below 1)', info['E_i_mean_left'] / phic, 0.05, 1.0)
ok &= check_range('7  maximum ion impact energy / (V0 + V0_lf)', emax / Vmax, 0.3, 1.1)

fig, ax = plt.subplots(2, 2, figsize=(13, 9))
a = ax[0, 0]
a.plot(t, sl * 1e2, label='left edge'); a.plot(t, L * 1e2 - sr * 1e2, label='right edge (distance from x = L)')
a2 = a.twinx(); a2.plot(t, expect, 'k:', lw=0.8); a2.set_ylabel('$\\phi(0,t)$ [V] (bin average)')
a.set_xlabel('t / T'); a.set_ylabel('sheath thickness [cm]'); a.legend(fontsize=8); a.set_title('sheath edges (n$_e$/n$_i$ = 0.5) over the RF period', fontsize=10)
a = ax[0, 1]
for p in range(0, nph, nph // 8):
    a.plot(x * 1e2, it_[p], lw=0.9, label='t/T = %.2f' % t[p])
a.set_xlabel('x [cm]'); a.set_ylabel('J$_e$ + J$_i$ + $\\epsilon_0\\,\\partial E/\\partial t$ [A/m$^2$]'); a.legend(fontsize=7, ncol=2)
a.set_title('total current is uniform in x (rms variation %.1f %% of the amplitude)' % (100 * spread / amp), fontsize=10)
a = ax[1, 0]
m = fp > 0
a.semilogy(e[m], fp[m], '.', ms=3, label='EEPF, central 10 %')
ee = np.linspace(0, 90, 200)
a.semilogy(ee, np.exp(np.polyval(np.polyfit(e[lo], np.log(fp[lo]), 1), ee)), 'C1--', label='1-8 eV slope: T = %.1f eV' % t_lo)
a.semilogy(ee, np.exp(np.polyval(np.polyfit(e[hi], np.log(fp[hi]), 1), ee)), 'C2--', label='25-80 eV slope: T = %.1f eV' % t_hi)
a.set_ylim(fp[m].min() * 0.5, fp[m].max() * 2); a.set_xlabel('electron energy [eV]'); a.set_ylabel('F/$\\sqrt{\\epsilon}$ [m$^{-3}$ eV$^{-3/2}$]'); a.legend(fontsize=8)
a.set_title('EEPF: two temperatures bracket 2/3 <$\\epsilon$> = %.1f eV' % teff, fontsize=10)
a = ax[1, 1]
a.plot(iedf[:, 0], iedf[:, 1], label='IEDF, x = 0')
a.axvline(info['E_i_mean_left'], color='gray', ls='--', label='mean %.0f eV' % info['E_i_mean_left'])
a.axvline(phic, color='r', ls=':', label='e <$\\phi_c$> = %.0f eV' % phic)
a.set_xlim(0, 1.05 * Vmax); a.set_xlabel('ion energy [eV]'); a.set_ylabel('m$^{-2}$ s$^{-1}$ eV$^{-1}$'); a.legend(fontsize=8)
a.set_title('IEDF: charge-exchange peaks, maximum %.0f eV' % emax, fontsize=10)
fig.suptitle('diagnostic self-consistency (T16): ' + os.path.abspath(d), fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(d, 'diag_checks.png'), dpi=120)
print('wrote', os.path.join(d, 'diag_checks.png'))
print('check_diag:', 'PASS' if ok else 'FAIL')
raise SystemExit(0 if ok else 1)
