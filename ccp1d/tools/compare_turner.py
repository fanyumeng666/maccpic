#!/usr/bin/env python3
"""compare_turner.py <run_dir> [more_run_dirs...]: time-averaged ion density of a
Case 1 run against the Turner et al. (2013) reference profile
(data/turner2013_case1_ni.dat, 129 nodes).  If the run grid differs from 129
nodes the reference is interpolated onto the run grid (as WarpX's
analysis_1d.py does).  Prints centre density, maximum and rms relative
deviation over all nodes and over the interior nodes, the left-right
asymmetry, and writes compare_turner.png in the first run directory."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import load_dat, load_info, profile_columns

here = os.path.dirname(os.path.abspath(__file__))
ref = np.loadtxt(os.path.join(here, '..', 'data', 'turner2013_case1_ni.dat'))
xr, nr = ref[:, 0], ref[:, 1]
dirs = sys.argv[1:] or ['.']
pc = profile_columns()

fig, ax = plt.subplots(2, 1, figsize=(9, 8), sharex=True, gridspec_kw={'height_ratios': [3, 1.4]})
ax[0].plot(xr * 1e2, nr, 'k-', lw=2.5, label='Turner et al. 2013, case 1 (5 codes agree to < 1 %)')
ok_all = True
for k, d in enumerate(dirs):
    prof, _ = load_dat(os.path.join(d, 'profiles.dat'))
    info = load_info(os.path.join(d, 'info.txt'))
    x, ni, ne = prof[:, pc['x']], prof[:, pc['n_i']], prof[:, pc['n_e']]
    nref = np.interp(x, xr, nr)
    mid = len(x) // 2
    nc = ni[mid - 2:mid + 3].mean()
    ncr = nr[len(nr) // 2 - 2:len(nr) // 2 + 3].mean()
    rel = (ni - nref) / nref
    inner = slice(1, len(x) - 1)
    asym = np.max(np.abs(ni - ni[::-1])) / nc
    nl, nrr = np.trapezoid(ni[:mid + 1], x[:mid + 1]), np.trapezoid(ni[mid:], x[mid:])
    label = os.path.basename(os.path.abspath(d)) or d
    print('==== %s  (nx = %d, nppc = %d, steps/period = %d, periods %d, seed %d)' %
          (label, info['nx'], info['nppc'], info['steps_per_period'], info['n_periods'], info['seed']))
    print('    centre n_i = %.4e  reference %.4e  deviation %+.2f %%' % (nc, ncr, 100 * (nc / ncr - 1)))
    print('    all nodes:      max |rel. dev.| = %.2f %%   rms = %.2f %%' % (100 * np.max(np.abs(rel)), 100 * np.sqrt(np.mean(rel**2))))
    print('    interior nodes: max |rel. dev.| = %.2f %%   rms = %.2f %%   mean = %+.2f %%' %
          (100 * np.max(np.abs(rel[inner])), 100 * np.sqrt(np.mean(rel[inner]**2)), 100 * np.mean(rel[inner])))
    print('    electrode nodes: run %.3e / %.3e, reference %.3e / %.3e' % (ni[0], ni[-1], nref[0], nref[-1]))
    print('    left-right: node-by-node max asymmetry %.2f %% of n_c, half-domain ion numbers differ by %.2f %%' %
          (100 * asym, 100 * abs(nl - nrr) / (nl + nrr)))
    print('    n_e centre / n_i centre = %.4f' % (ne[mid - 2:mid + 3].mean() / nc))
    ok = abs(nc / ncr - 1) < 0.05 and np.max(np.abs(rel[inner])) < 0.10 and np.sqrt(np.mean(rel[inner]**2)) < 0.05
    print('    T14 criteria (centre < 5 %%, interior max < 10 %%, interior rms < 5 %%): %s' % ('PASS' if ok else 'FAIL'))
    ok_all &= ok
    ax[0].plot(x * 1e2, ni, '-', lw=1.2, label='%s: n_i (centre %+.1f %%)' % (label, 100 * (nc / ncr - 1)))
    if k == 0:
        ax[0].plot(x * 1e2, ne, '--', lw=0.9, color='C0', alpha=0.7, label=label + ': n_e')
    ax[1].plot(x * 1e2, 100 * rel, '-', lw=1)
ax[0].set_ylabel('density [m$^{-3}$]'); ax[0].legend(fontsize=8); ax[0].set_title('Turner 2013 case 1: He, 450 V, 13.56 MHz, 6.7 cm, 30 mTorr; time average of the last %d periods' % int(info['n_avg']), fontsize=10)
ax[0].set_ylim(bottom=0)
ax[1].axhline(0, color='k', lw=0.8); ax[1].axhspan(-5, 5, color='gray', alpha=0.15, label='$\\pm$5 %')
ax[1].set_xlabel('x [cm]'); ax[1].set_ylabel('(n$_i$ - ref)/ref [%]'); ax[1].set_ylim(-25, 25); ax[1].legend(fontsize=8)
fig.tight_layout()
out = os.path.join(dirs[0], 'compare_turner.png')
fig.savefig(out, dpi=130)
print('wrote', out)
raise SystemExit(0 if ok_all else 1)
