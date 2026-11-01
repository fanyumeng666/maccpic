#!/usr/bin/env python3
"""plot_profiles.py <run_dir>: time-averaged profiles from profiles.dat
(densities, potential, field, electron mean energy and T_x, ionization rate,
J.E heating of electrons and ions, ion flux) -> profiles.png."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import load_dat, load_info, profile_columns, run_dir_from_argv

d = run_dir_from_argv()
prof, _ = load_dat(os.path.join(d, 'profiles.dat'))
info = load_info(os.path.join(d, 'info.txt'))
pc = profile_columns()
x = prof[:, pc['x']] * 1e2
fig, ax = plt.subplots(2, 4, figsize=(18, 8))
a = ax[0, 0]; a.plot(x, prof[:, pc['n_e']], label='n$_e$'); a.plot(x, prof[:, pc['n_i']], label='n$_i$')
a.set_ylabel('m$^{-3}$'); a.legend(); a.set_title('densities'); a.set_ylim(bottom=0)
a = ax[0, 1]; a.plot(x, prof[:, pc['phi']]); a.set_ylabel('V'); a.set_title('time-averaged potential')
a = ax[0, 2]; a.plot(x, prof[:, pc['E']]); a.set_ylabel('V/m'); a.set_title('time-averaged field')
a = ax[0, 3]; a.plot(x, prof[:, pc['eps_e']], label=r'<$\epsilon_e$>'); a.plot(x, 2 * prof[:, pc['eps_x_e']], '--', label=r'T$_{x,e}$ = 2<$\epsilon_{x,e}$>')
a.plot(x, 2.0 / 3.0 * prof[:, pc['eps_e']], ':', label=r'T$_{eff}$ = 2/3 <$\epsilon_e$>')
a.set_ylabel('eV'); a.legend(fontsize=8); a.set_title('electron energy'); a.set_ylim(bottom=0)
a = ax[1, 0]; a.plot(x, prof[:, pc['ioniz']]); a.set_ylabel('m$^{-3}$ s$^{-1}$'); a.set_title('ionization rate'); a.set_ylim(bottom=0)
a = ax[1, 1]; a.plot(x, prof[:, pc['P_work_e']], label='electrons (exact work deposit)'); a.plot(x, prof[:, pc['P_work_i']], label='ions')
a.plot(x, prof[:, pc['PeE']], ':', color='C0', label='<J$_e$.E> grid product'); a.plot(x, prof[:, pc['PiE']], ':', color='C1', label='<J$_i$.E>')
a.axhline(0, color='k', lw=0.6); a.set_ylabel('W/m$^3$'); a.legend(fontsize=7); a.set_title('time-averaged heating power density')
a = ax[1, 2]; a.plot(x, prof[:, pc['flux_i']]); a.axhline(0, color='k', lw=0.6); a.set_ylabel('m$^{-2}$ s$^{-1}$'); a.set_title('ion flux n$_i$ v$_i$')
a = ax[1, 3]; a.plot(x, prof[:, pc['eps_i']]); a.set_ylabel('eV'); a.set_title('mean ion energy'); a.set_yscale('log')
for row in ax:
    for a in row:
        a.set_xlabel('x [cm]')
fig.suptitle('%s: time averages over the last %d periods (nx = %d, nppc = %d)' % (os.path.abspath(d), int(info['n_avg']), int(info['nx']), int(info['nppc'])), fontsize=10)
fig.tight_layout()
out = os.path.join(d, 'profiles.png')
fig.savefig(out, dpi=120)
print('wrote', out)
