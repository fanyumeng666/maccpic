#!/usr/bin/env python3
"""plot_dist.py <run_dir>: ion energy distributions at both electrodes (IEDF),
the electron flux energy distribution at the electrodes (EEDF) and the bulk
EEPF on a semi-log scale with the Maxwellian slope -> dist.png."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import load_dat, load_info, run_dir_from_argv

d = run_dir_from_argv()
info = load_info(os.path.join(d, 'info.txt'))
il, _ = load_dat(os.path.join(d, 'iedf_left.dat'))
ir, _ = load_dat(os.path.join(d, 'iedf_right.dat'))
el, _ = load_dat(os.path.join(d, 'eedf_left.dat'))
ep, _ = load_dat(os.path.join(d, 'eepf.dat'))
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
a = ax[0]
a.plot(il[:, 0], il[:, 1], label='x = 0 (driven)')
a.plot(ir[:, 0], ir[:, 1], label='x = L (grounded)')
mean_l = np.sum(il[:, 0] * il[:, 1]) / max(np.sum(il[:, 1]), 1e-300)
a.axvline(mean_l, color='gray', ls='--', lw=0.8, label='mean %.0f eV' % mean_l)
a.set_xlim(0, 1.1 * info['V0']); a.set_xlabel('ion energy [eV]'); a.set_ylabel('flux per eV [m$^{-2}$ s$^{-1}$ eV$^{-1}$]')
a.set_title('IEDF at the electrodes (V$_0$ = %.0f V)' % info['V0']); a.legend(fontsize=8)
a = ax[1]
m = el[:, 1] > 0
a.semilogy(el[m, 0], el[m, 1], '.', ms=3, label='x = 0')
er, _ = load_dat(os.path.join(d, 'eedf_right.dat'))
m = er[:, 1] > 0
a.semilogy(er[m, 0], er[m, 1], '.', ms=3, label='x = L')
a.set_xlim(0, 100); a.set_xlabel('electron energy [eV]'); a.set_ylabel('flux per eV'); a.set_title('electron flux energy distribution at the electrodes'); a.legend(fontsize=8)
a = ax[2]
m = ep[:, 2] > 0
a.semilogy(ep[m, 0], ep[m, 2], '.', ms=3, label='central 10 %')
# Maxwellian with the centre T_eff for reference
te = info.get('T_eff_center', 3.0)
ne = info.get('n_e_center', 1e14)
e = ep[:, 0]
a.semilogy(e, ne * 2 / np.sqrt(np.pi) * te**-1.5 * np.exp(-e / te), 'k--', lw=1, label='Maxwellian, T = %.2f eV' % te)
a.set_xlim(0, info['eepf_emax']); a.set_ylim(bottom=max(ep[m, 2].min() * 0.5, 1e-300) if np.any(m) else 1)
a.set_xlabel('electron energy [eV]'); a.set_ylabel('f$_p$ = F/$\\sqrt{\\epsilon}$ [m$^{-3}$ eV$^{-3/2}$]'); a.set_title('EEPF (Maxwellian = straight line)'); a.legend(fontsize=8)
fig.suptitle(os.path.abspath(d), fontsize=9)
fig.tight_layout()
out = os.path.join(d, 'dist.png')
fig.savefig(out, dpi=120)
print('wrote', out)
