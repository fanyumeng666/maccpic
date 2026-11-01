#!/usr/bin/env python3
"""check_sheath.py <run_dir>: T11, the collisionless sheath test (inputs/sheath.in).

Hard checks: finite data, monotone particle loss, the exact bookkeeping identity,
total energy conservation, the discrete Gauss law, left-right symmetry, flux
balance Gamma_e = Gamma_i per electrode, constancy of the perpendicular electron
temperature (nothing may touch vy, vz without collisions).  Loose checks: the
plasma potential against the floating potential built with the x-temperature
(in a collisionless 1D plasma T_x collapses while T_perp stays, so the usual
3.5 T_e does not apply), and the IEDF against e phi_c.  Writes sheath.png."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import (load_dat, load_info, history_columns, profile_columns, run_dir_from_argv,
                    check, check_range, QE, ME, MION, KB)

d = run_dir_from_argv()
h, _ = load_dat(os.path.join(d, 'history.dat'))
c = history_columns()
prof, _ = load_dat(os.path.join(d, 'profiles.dat'))
pc = profile_columns()
info = load_info(os.path.join(d, 'info.txt'))
ok = True

# (i) finite, particle numbers never increase
ok &= check('T11.i non-finite entries in history', float(np.sum(~np.isfinite(h))), 0)
ok &= check('T11.i increases of N_e+N_i between periods', float(np.sum(np.diff(h[:, c['N_e']] + h[:, c['N_i']]) > 0)), 0)
# (ii) energy: exact identity and physical conservation
ke = h[:, c['KE_e']] + h[:, c['KE_i']]
ident = ke - ke[0] - (h[:, c['Work_cum']] - h[:, c['E_wall_cum']])
ok &= check('T11.ii exact identity |KE-KE0-(Work-Ewall)|/KE0', np.max(np.abs(ident)) / ke[0], 1e-6)   # round-off over ~1e10 pushes
tot = ke + h[:, c['W_field']] + h[:, c['E_wall_cum']]
drift = (tot.max() - tot.min()) / ke[0]
ok &= check('T11.ii total energy spread / KE(0)', drift, 1e-2)
print('        energy change over the run: %+.3e of KE(0) (%s); %.0f %% of the initial energy went to the walls' %
      ((tot[-1] - tot[0]) / ke[0], 'heating' if tot[-1] > tot[0] else 'cooling', 100 * h[-1, c['E_wall_cum']] / ke[0]))
# (iii) Gauss law
ok &= check('T11.iii max Gauss residual', np.max(h[1:, c['gauss_res_max']]), 1e-10)
# (iv) symmetry of the window-averaged ion density and potential
x = prof[:, pc['x']]
ni = prof[:, pc['n_i']]
phi = prof[:, pc['phi']]
mid = len(x) // 2
nc = ni[mid - 2:mid + 3].mean()
# node-by-node the window average still carries the ion shot noise (ions hardly move
# within the window), so compare the two halves integrated: N_left vs N_right
nl, nr = np.trapezoid(ni[:mid + 1], x[:mid + 1]), np.trapezoid(ni[mid:], x[mid:])
ok &= check('T11.iv |N_i,left - N_i,right| / (N_i,left + N_i,right)', abs(nl - nr) / (nl + nr), 0.03)
phic = phi[mid - 2:mid + 3].mean()
pl, pr = phi[1:mid].mean(), phi[mid + 1:-1].mean()
ok &= check('T11.iv |<phi>_left/<phi>_right - 1|', abs(pl / pr - 1), 0.1)
print('        node-by-node max asymmetry: n_i %.1f %% of n_c (shot noise per node %.1f %%), phi %.1f %% of phi_c' %
      (100 * np.max(np.abs(ni - ni[::-1])) / nc, 100 / np.sqrt(nc * info['dx'] / info['weight']), 100 * np.max(np.abs(phi - phi[::-1])) / phic))
# (v) perpendicular electron temperature must not change (vy, vz untouched)
tp0, tp1 = h[0, c['Tperp_e']], h[-1, c['Tperp_e']]
ok &= check('T11.v |Tperp_e(end)/Tperp_e(0) - 1|', abs(tp1 / tp0 - 1), 0.02)
tx0, tx1 = h[0, c['Tx_e']], h[-1, c['Tx_e']]
print('        electrons: T_x %.3f -> %.3f eV, T_perp %.3f -> %.3f eV (collisionless 1D: T_x collapses, T_perp conserved)' % (tx0, tx1, tp0, tp1))
# (vi) flux balance per electrode in the window
for side in ['left', 'right']:
    gi, ge = info['Gamma_i_' + side], info['Gamma_e_' + side]
    ok &= check('T11.vi |Gamma_e/Gamma_i - 1| ' + side, abs(ge / gi - 1), 0.1)
ok &= check('T11.vi |Gamma_i_left/Gamma_i_right - 1|', abs(info['Gamma_i_left'] / info['Gamma_i_right'] - 1), 0.1)
# (vii) plasma potential vs floating potential with T_x (loose), IEDF vs e phi_c (loose)
txc = info.get('Tx_e_center', tx1)
vf = 0.5 * txc * np.log(MION / (2 * np.pi * ME))
print('        centre: T_x,e = %.3f eV, phi_c = %.3f V, (T_x/2) ln(M/2 pi m_e) = %.3f V, phi_c/T_x = %.2f' % (txc, phic, vf, phic / txc))
ok &= check_range('T11.vii phi_c / V_float(T_x)', phic / vf, 0.2, 2.0)
iedf, _ = load_dat(os.path.join(d, 'iedf_left.dat'))
e_bin, f = iedf[:, 0], iedf[:, 1]
emean = np.sum(e_bin * f) / np.sum(f)
print('        IEDF left: peak at %.2f eV, mean %.2f eV (ions were accelerated earlier, when phi_c was higher)' % (e_bin[np.argmax(f)], emean))
ok &= check_range('T11.vii IEDF mean energy / e phi_c', emean / phic, 0.3, 5.0)
cs = np.sqrt(txc * QE / MION)
print('        Bohm flux n_c c_s(T_x) = %.3e, measured Gamma_i = %.3e m^-2 s^-1' % (nc * cs, info['Gamma_i_left']))

fig, ax = plt.subplots(2, 3, figsize=(16, 8.5))
a = ax[0, 0]
a.plot(x * 1e2, prof[:, pc['n_e']], label='n$_e$')
a.plot(x * 1e2, ni, label='n$_i$')
a.set_xlabel('x [cm]'); a.set_ylabel('m$^{-3}$'); a.legend(); a.set_title('density, window average (last %d periods)' % int(info['n_avg']))
a = ax[0, 1]
a.plot(x * 1e2, phi, label=r'$\phi$')
a.axhline(vf, color='r', ls='--', label=r'$(T_x/2)\ln(M/2\pi m_e)$ = %.2f V' % vf)
a.set_xlabel('x [cm]'); a.set_ylabel(r'$\phi$ [V]'); a.legend(); a.set_title('potential (walls at 0 V)')
a = ax[0, 2]
t = h[:, c['t']] * 1e6
a.plot(t, h[:, c['Tx_e']], label='T$_{x,e}$')
a.plot(t, h[:, c['Tperp_e']], label=r'T$_{\perp,e}$')
a.plot(t, h[:, c['Tx_i']], label='T$_{x,i}$')
a.plot(t, h[:, c['Ei_L']], 'k--', label='mean ion impact energy (left)')
a.set_xlabel('t [$\\mu$s]'); a.set_ylabel('eV'); a.legend(fontsize=8); a.set_title('temperatures: T$_x$ collapses without collisions')
a = ax[1, 0]
a.plot(t, h[:, c['KE_e']], label='KE$_e$')
a.plot(t, h[:, c['KE_i']], label='KE$_i$')
a.plot(t, h[:, c['W_field']], label='W$_{field}$')
a.plot(t, h[:, c['E_wall_cum']], label='to walls (cum.)')
a.plot(t, tot, 'k', label='total')
a.set_xlabel('t [$\\mu$s]'); a.set_ylabel('J m$^{-2}$'); a.legend(fontsize=8); a.set_title('energy budget: total conserved to %.1e' % drift)
a = ax[1, 1]
a.plot(t, h[:, c['N_e']], label='N$_e$')
a.plot(t, h[:, c['N_i']], label='N$_i$')
a.set_xlabel('t [$\\mu$s]'); a.set_ylabel('macro-particles'); a.legend(); a.set_title('particle loss (ion transit L/2c$_s$ = %.1f $\\mu$s)' % (info['ion_transit_half_gap'] * 1e6))
a = ax[1, 2]
a.plot(e_bin, f, label='IEDF left')
iedr, _ = load_dat(os.path.join(d, 'iedf_right.dat'))
a.plot(iedr[:, 0], iedr[:, 1], label='IEDF right')
a.axvline(phic, color='r', ls='--', label=r'e$\phi_c$')
a.set_xlim(0, 5); a.set_xlabel('energy [eV]'); a.set_ylabel('m$^{-2}$ s$^{-1}$ eV$^{-1}$'); a.legend(); a.set_title('ion energy distribution at the walls')
fig.suptitle('collisionless sheath test (T11): ' + os.path.abspath(d), fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(d, 'sheath.png'), dpi=120)
print('wrote', os.path.join(d, 'sheath.png'))
print('check_sheath:', 'PASS' if ok else 'FAIL')
raise SystemExit(0 if ok else 1)
