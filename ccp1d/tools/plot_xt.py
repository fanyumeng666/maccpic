#!/usr/bin/env python3
"""plot_xt.py <run_dir>: phase-resolved (x, t/T) maps from xt_*.dat:
potential, electron and ion density, electron and ion current density,
electron heating (exact work deposit and the grid product J_e.E) and the
ionization rate.  The sheath edges (n_e/n_i = 0.5) are overlaid in white.
Writes xt.png."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import load_dat, load_info, run_dir_from_argv

d = run_dir_from_argv()
info = load_info(os.path.join(d, 'info.txt'))
L, V0 = info['L'], info['V0']


def xt(name):
    a, _ = load_dat(os.path.join(d, 'xt_%s.dat' % name))
    return a


phi, ne, ni = xt('phi'), xt('ne'), xt('ni')
je, ji, pee, ioniz, pwork = xt('Je'), xt('Ji'), xt('PeE'), xt('ioniz'), xt('Pe_work')
nph, nx1 = phi.shape
x = np.linspace(0, L, nx1) * 1e2
t = (np.arange(nph) + 0.5) / nph


def sheath_edges(ne, ni):
    r = ne / np.maximum(ni, 1e-300)
    left = np.zeros(nph); right = np.zeros(nph)
    for p in range(nph):
        above = np.nonzero(r[p] >= 0.5)[0]
        left[p] = x[above[0]] if len(above) else x[-1] / 2
        right[p] = x[above[-1]] if len(above) else x[-1] / 2
    return left, right


sl, sr = sheath_edges(ne, ni)
fig, ax = plt.subplots(2, 4, figsize=(20, 9))
panels = [
    (ax[0, 0], phi, 'potential $\\phi$ [V]', 'viridis', False),
    (ax[0, 1], ne, 'n$_e$ [m$^{-3}$]', 'viridis', False),
    (ax[0, 2], ni, 'n$_i$ [m$^{-3}$]', 'viridis', False),
    (ax[0, 3], je, 'J$_e$ [A/m$^2$]', 'RdBu_r', True),
    (ax[1, 0], ji, 'J$_i$ [A/m$^2$]', 'RdBu_r', True),
    (ax[1, 1], pwork, 'electron heating, exact work deposit [W/m$^3$]', 'RdBu_r', True),
    (ax[1, 2], pee, 'grid product J$_e$.E [W/m$^3$]', 'RdBu_r', True),
    (ax[1, 3], ioniz, 'ionization rate [m$^{-3}$ s$^{-1}$]', 'magma', False),
]
for a, z, title, cmap, sym in panels:
    if sym:
        vmax = np.percentile(np.abs(z), 99.5)
        im = a.pcolormesh(x, t, z, cmap=cmap, vmin=-vmax, vmax=vmax, shading='auto')
    else:
        im = a.pcolormesh(x, t, z, cmap=cmap, shading='auto')
    fig.colorbar(im, ax=a, pad=0.02)
    a.plot(sl, t, 'w--', lw=0.8); a.plot(sr, t, 'w--', lw=0.8)
    a.set_title(title, fontsize=10); a.set_xlabel('x [cm]'); a.set_ylabel('t / T')
ax[0, 0].text(0.02, 0.97, 'white dashed: sheath edge (n$_e$/n$_i$ = 0.5)', transform=ax[0, 0].transAxes, color='w', fontsize=8, va='top')
fig.suptitle('%s: phase-resolved maps, average over the last %d periods; $\\phi(0,t) = %.0f\\,\\sin(2\\pi t/T)$ V' % (os.path.abspath(d), int(info['n_avg']), V0), fontsize=10)
fig.tight_layout()
out = os.path.join(d, 'xt.png')
fig.savefig(out, dpi=110)
print('wrote', out)
print('sheath edge excursion: left %.2f-%.2f cm (widest at t/T = %.2f), right %.2f-%.2f cm (widest at t/T = %.2f)' %
      (sl.min(), sl.max(), t[np.argmax(sl)], sr.min(), sr.max(), t[np.argmin(sr)]))
