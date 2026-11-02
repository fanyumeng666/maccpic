#!/usr/bin/env python3
"""dual_frequency.py <run_dir_Vlf0> <run_dir_1> [<run_dir_2> ...]: T17, the dual-frequency
extension.  For each run: the low-frequency amplitude V_lf, the ion energy distribution at
the driven electrode (x = 0), its mean, its 10 %-90 % width and its maximum, the centre
density and the absorbed power.  Criterion: the IEDF width grows monotonically with V_lf
and the maximum ion energy follows V0 + V_lf.  Writes dual_frequency.png in the first dir."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import load_dat, load_info, profile_columns

dirs = sys.argv[1:]
if len(dirs) < 2:
    sys.exit('usage: dual_frequency.py <dir_Vlf0> <dir1> [...]')
pc = profile_columns()
rows = []
for d in dirs:
    info = load_info(os.path.join(d, 'info.txt'))
    iedf, _ = load_dat(os.path.join(d, 'iedf_left.dat'))
    e, f = iedf[:, 0], iedf[:, 1]
    de = e[1] - e[0]
    flux = np.sum(f) * de
    cdf = np.cumsum(f) * de / flux
    e10, e90 = np.interp(0.1, cdf, e), np.interp(0.9, cdf, e)
    emean = np.sum(e * f) * de / flux
    emax = e[np.nonzero(f)[0][-1]] + 0.5 * de
    prof, _ = load_dat(os.path.join(d, 'profiles.dat'))
    mid = len(prof) // 2
    vlf = info.get('V0_lf', 0.0) if info['freq'] > info.get('freq_lf', 0.0) else info['V0']
    vhf = info['V0'] if info['freq'] > info.get('freq_lf', 0.0) else info.get('V0_lf', 0.0)
    nions = flux * info['window_duration'] / info['weight']          # macro-ions in the histogram
    rows.append(dict(name=os.path.basename(os.path.abspath(d)) or d, vlf=vlf, vhf=vhf, e=e, f=f, flux=flux, emean=emean,
                     width=e90 - e10, e10=e10, e90=e90, emax=emax, nc=prof[mid - 2:mid + 3, pc['n_i']].mean(),
                     sig=1.5 * (e90 - e10) / np.sqrt(max(nions, 1.0)),   # statistical noise of a quantile width
                     pvi=info['P_VI'], gi=info['Gamma_i_left'], te=info['T_eff_center'], phic=prof[mid - 2:mid + 3, pc['phi']].mean()))
rows.sort(key=lambda r: r['vlf'])
print('%-14s %6s %8s %8s %8s %12s %9s %11s %9s %7s %8s' % ('run', 'V_lf', '<E_i>', 'E10', 'E90', 'width(+-)', 'E_max', 'n_i centre', 'Gamma_i', 'P_VI', '<phi_c>'))
for r in rows:
    print('%-14s %6.0f %8.1f %8.1f %8.1f %7.1f(%4.1f) %9.0f %11.3e %9.3e %7.1f %8.1f' %
          (r['name'], r['vlf'], r['emean'], r['e10'], r['e90'], r['width'], r['sig'], r['emax'], r['nc'], r['gi'], r['pvi'], r['phic']))
widths = [r['width'] for r in rows]; sigs = [r['sig'] for r in rows]
# monotonic within the noise: no decrease larger than 2 sigma, and a net increase larger than 3 sigma
no_drop = all(widths[i + 1] >= widths[i] - 2 * np.hypot(sigs[i], sigs[i + 1]) for i in range(len(widths) - 1))
net = widths[-1] - widths[0] > 3 * np.hypot(sigs[0], sigs[-1])
ok = no_drop and net
print('T17 IEDF 10-90 %% width grows with V_lf (no step down beyond 2 sigma, net rise beyond 3 sigma): %s  (net %.1f eV, sigma %.1f eV)' %
      ('PASS' if ok else 'FAIL', widths[-1] - widths[0], np.hypot(sigs[0], sigs[-1])))
emaxs = [r['emax'] for r in rows]
print('    E_max monotonic: %s (%s)' % ('yes' if all(emaxs[i + 1] > emaxs[i] for i in range(len(emaxs) - 1)) else 'no', ', '.join('%.0f' % q for q in emaxs)))
ratios = [r['emax'] / (r['vhf'] + r['vlf']) for r in rows]
print('    E_max / (V0 + V_lf): ' + ', '.join('%.2f' % q for q in ratios))
ok2 = all(0.3 <= q <= 1.1 for q in ratios)
print('T17 E_max within 0.3..1.1 of V0 + V_lf: %s' % ('PASS' if ok2 else 'FAIL'))
base = rows[0]
print('    centre density change relative to V_lf = 0: ' + ', '.join('%+.1f %%' % (100 * (r['nc'] / base['nc'] - 1)) for r in rows[1:]))

fig, ax = plt.subplots(1, 3, figsize=(17, 5))
for r in rows:
    ax[0].plot(r['e'], r['f'] / r['flux'], lw=1, label='V$_{lf}$ = %.0f V: mean %.0f eV, 10-90 %% width %.0f eV' % (r['vlf'], r['emean'], r['width']))
ax[0].set_xlim(0, 1.05 * max(r['emax'] for r in rows)); ax[0].set_xlabel('ion energy [eV]'); ax[0].set_ylabel('IEDF / flux [eV$^{-1}$]')
ax[0].set_title('IEDF at the driven electrode, 13.56 MHz 450 V + 1.356 MHz V$_{lf}$', fontsize=10); ax[0].legend(fontsize=8)
v = [r['vlf'] for r in rows]
ax[1].plot(v, [r['width'] for r in rows], 'o-', label='10-90 % width')
ax[1].plot(v, [r['emean'] for r in rows], 's-', label='mean energy')
ax[1].plot(v, [r['emax'] for r in rows], '^-', label='maximum energy')
ax[1].plot(v, [r['vhf'] + r['vlf'] for r in rows], 'k:', label='V$_0$ + V$_{lf}$')
ax[1].set_xlabel('V$_{lf}$ [V]'); ax[1].set_ylabel('eV'); ax[1].legend(fontsize=8); ax[1].set_title('ion energies versus the low-frequency amplitude', fontsize=10)
ax[2].plot(v, [r['nc'] for r in rows], 'o-', label='n$_i$ centre')
ax2 = ax[2].twinx(); ax2.plot(v, [r['pvi'] for r in rows], 's--', color='C3', label='absorbed power'); ax2.set_ylabel('W/m$^2$')
ax[2].set_xlabel('V$_{lf}$ [V]'); ax[2].set_ylabel('m$^{-3}$'); ax[2].set_title('density and power: ion energy is controlled separately', fontsize=10)
h1, l1 = ax[2].get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax[2].legend(h1 + h2, l1 + l2, fontsize=8)
fig.tight_layout()
out = os.path.join(dirs[0], 'dual_frequency.png')
fig.savefig(out, dpi=130)
print('wrote', out)
raise SystemExit(0 if (ok and ok2) else 1)
