#!/usr/bin/env python3
"""plot_history.py <run_dir>: particle numbers, densities, fluxes, energies and
power versus period from history.dat -> history.png; prints the steady-state
indicators (drift of the particle numbers over the last 20 % of the run)."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import load_dat, history_columns, run_dir_from_argv

d = run_dir_from_argv()
h, _ = load_dat(os.path.join(d, 'history.dat'))
c = history_columns()
p = h[:, c['period']]
t = h[:, c['t']] * 1e6

fig, ax = plt.subplots(2, 3, figsize=(14, 8))
a = ax[0, 0]
a.plot(p, h[:, c['N_e']], label='N_e')
a.plot(p, h[:, c['N_i']], label='N_i')
a.set_xlabel('period'); a.set_ylabel('macro-particles (all ranks)'); a.legend(); a.set_title('particle numbers')
a = ax[0, 1]
a.plot(p, h[:, c['n_e_center']], label='n_e centre')
a.plot(p, h[:, c['n_i_center']], label='n_i centre')
a.plot(p, h[:, c['n_e_mean']], '--', label='<n_e>')
a.plot(p, h[:, c['n_i_mean']], '--', label='<n_i>')
a.set_xlabel('period'); a.set_ylabel('m$^{-3}$'); a.legend(); a.set_title('densities')
a = ax[0, 2]
for k in ['Gi_L', 'Gi_R', 'Ge_L', 'Ge_R']:
    a.plot(p, h[:, c[k]], label=k)
a.set_xlabel('period'); a.set_ylabel('m$^{-2}$ s$^{-1}$'); a.legend(); a.set_title('wall fluxes (period average)')
a = ax[1, 0]
ke = h[:, c['KE_e']] + h[:, c['KE_i']]
tot = ke + h[:, c['W_field']] + h[:, c['E_wall_cum']]
a.plot(p, h[:, c['KE_e']], label='KE_e')
a.plot(p, h[:, c['KE_i']], label='KE_i')
a.plot(p, h[:, c['W_field']], label='W_field')
a.plot(p, h[:, c['E_wall_cum']], label='to walls (cum.)')
a.plot(p, tot, 'k', label='total')
a.set_xlabel('period'); a.set_ylabel('J m$^{-2}$'); a.legend(fontsize=8); a.set_title('energies')
a = ax[1, 1]
a.plot(p, h[:, c['P_e']], label='P_e (work on electrons)')
a.plot(p, h[:, c['P_i']], label='P_i (work on ions)')
a.plot(p, h[:, c['P_VI']], 'k--', label='P_VI (source)')
a.set_xlabel('period'); a.set_ylabel('W m$^{-2}$'); a.legend(fontsize=8); a.set_title('power (period average)')
a = ax[1, 2]
a.semilogy(p[1:], np.maximum(h[1:, c['gauss_res_max']], 1e-18), label='Gauss residual (max)')
if h[-1, c['t']] > 0:
    rel = np.abs(tot - tot[0]) / max(abs(tot[0]), 1e-300)
    a.semilogy(p[1:], np.maximum(rel[1:], 1e-18), label='|E_tot - E_tot(0)|/E_tot(0)')
a.set_xlabel('period'); a.legend(fontsize=8); a.set_title('conservation checks')
fig.suptitle(os.path.abspath(d))
fig.tight_layout()
out = os.path.join(d, 'history.png')
fig.savefig(out, dpi=120)
print('wrote', out)

n = len(p)
if n > 10:
    k = max(2, n // 5)
    for name in ['N_e', 'N_i', 'n_i_center']:
        y = h[-k:, c[name]]
        slope = np.polyfit(p[-k:], y, 1)[0]
        print('%-12s last 20%%: mean %.4e, drift %.2f %% per 100 periods' % (name, y.mean(), 100 * slope * 100 / y.mean()))
print('final: N_e %.0f  N_i %.0f  n_i_center %.4e  wall %.1f s' % (h[-1, c['N_e']], h[-1, c['N_i']], h[-1, c['n_i_center']], h[-1, c['wall']]))
