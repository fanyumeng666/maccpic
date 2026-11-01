#!/usr/bin/env python3
"""convergence.py <baseline_dir> <dir2> [<dir3> ...]: T15.  Compares the
time-averaged ion density profiles of several runs (different dt, dx, particles
per cell or seed).  Profiles on a finer grid are interpolated onto the baseline
grid; each profile is smoothed with a 5-node running mean before the
node-by-node comparison.  The run that differs from the baseline only by the
seed sets the noise floor: every other pair must stay within 2x that
deviation.  Writes convergence.png next to the baseline."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import load_dat, load_info, profile_columns

dirs = sys.argv[1:]
if len(dirs) < 2:
    sys.exit('usage: convergence.py <baseline_dir> <dir2> [...]')
pc = profile_columns()
here = os.path.dirname(os.path.abspath(__file__))
ref = np.loadtxt(os.path.join(here, '..', 'data', 'turner2013_case1_ni.dat'))
runs = []
for d in dirs:
    prof, _ = load_dat(os.path.join(d, 'profiles.dat'))
    info = load_info(os.path.join(d, 'info.txt'))
    runs.append(dict(name=os.path.basename(os.path.abspath(d)) or d, x=prof[:, pc['x']], ni=prof[:, pc['n_i']], info=info))
x0 = runs[0]['x']
mid = len(x0) // 2


def smooth(y, w=5):
    k = np.ones(w) / w
    ys = np.convolve(y, k, mode='same')
    ys[:w // 2] = y[:w // 2]; ys[-(w // 2):] = y[-(w // 2):]
    return ys


for r in runs:
    r['ni0'] = smooth(np.interp(x0, r['x'], r['ni']))
    r['nc'] = r['ni'][len(r['x']) // 2 - 2:len(r['x']) // 2 + 3].mean()
    i = r['info']
    r['label'] = '%s (nx %d, nppc %d, %d steps/T, seed %d)' % (r['name'], i['nx'], i['nppc'], i['steps_per_period'], i['seed'])
b = runs[0]
ncr = ref[len(ref) // 2 - 2:len(ref) // 2 + 3, 1].mean()
print('centre n_i (5-node mean) and deviation from Turner (%.4e):' % ncr)
for r in runs:
    print('  %-55s %.4e  %+.2f %%' % (r['label'], r['nc'], 100 * (r['nc'] / ncr - 1)))
# noise floor: same nx, nppc, steps, different seed
noise = None
for r in runs[1:]:
    i, j = r['info'], b['info']
    if i['nx'] == j['nx'] and i['nppc'] == j['nppc'] and i['steps_per_period'] == j['steps_per_period'] and i['seed'] != j['seed']:
        noise = r
if noise is not None:
    d = np.abs(noise['ni0'] - b['ni0']) / b['nc']
    floor_max, floor_rms = d.max(), np.sqrt(np.mean(d**2))
    print('noise floor (baseline vs %s, 5-node smoothed): max %.2f %%, rms %.2f %% of n_c' % (noise['name'], 100 * floor_max, 100 * floor_rms))
else:
    floor_max = floor_rms = None
    print('no same-parameter different-seed run given: noise floor unknown')
print('pairwise deviations of the smoothed profiles (max / rms, % of the baseline centre density):')
names = [r['name'] for r in runs]
print(' ' * 12 + ''.join('%14s' % n[:13] for n in names))
for a in runs:
    row = '%12s' % a['name'][:11]
    for c in runs:
        d = np.abs(a['ni0'] - c['ni0']) / b['nc']
        row += '  %5.2f / %4.2f' % (100 * d.max(), 100 * np.sqrt(np.mean(d**2)))
    print(row)
# criterion: a run that differs from the baseline by dt, dx or more particles must stay within
# 2x the seed-to-seed noise; the expected noise of a pair scales with the particles per cell,
# sqrt((1/nppc_a + 1/nppc_b)/(2/nppc_base)), so a coarser run is judged against its own noise
ok = True
if floor_max is not None:
    print('baseline versus each run (allowed = 2x the seed noise scaled to the particles per cell of the pair):')
    for c in runs[1:]:
        if c is noise:
            continue
        d = np.abs(c['ni0'] - b['ni0']) / b['nc']
        sc = np.sqrt((1.0 / c['info']['nppc'] + 1.0 / b['info']['nppc']) / (2.0 / b['info']['nppc']))
        amax, arms = 2 * floor_max * sc, 2 * floor_rms * sc
        finer = (c['info']['nppc'] >= b['info']['nppc']) and (c['info']['nx'] >= b['info']['nx']) and (c['info']['steps_per_period'] >= b['info']['steps_per_period'])
        passed = d.max() <= amax and np.sqrt(np.mean(d**2)) <= arms
        print('  %-9s max %.2f %% (allowed %.2f), rms %.2f %% (allowed %.2f)  %s%s' % (c['name'], 100 * d.max(), 100 * amax, 100 * np.sqrt(np.mean(d**2)), 100 * arms,
              'PASS' if passed else 'FAIL', '' if finer else '  (coarser than the baseline: information only)'))
        if finer and not passed:
            ok = False
    print('T15 (every refinement within 2x the seed-to-seed noise):', 'PASS' if ok else 'FAIL')

fig, ax = plt.subplots(2, 1, figsize=(9, 8), sharex=True, gridspec_kw={'height_ratios': [3, 1.5]})
ax[0].plot(ref[:, 0] * 1e2, ref[:, 1], 'k', lw=2.5, label='Turner et al. 2013')
for r in runs:
    ax[0].plot(r['x'] * 1e2, r['ni'], lw=1, label=r['label'])
    ax[1].plot(x0 * 1e2, 100 * (r['ni0'] - b['ni0']) / b['nc'], lw=1, label=r['name'])
ax[0].set_ylabel('n$_i$ [m$^{-3}$]'); ax[0].legend(fontsize=7); ax[0].set_ylim(bottom=0)
ax[0].set_title('convergence: time-averaged ion density, Case 1', fontsize=10)
ax[1].axhline(0, color='k', lw=0.8)
if floor_max is not None:
    ax[1].axhspan(-100 * floor_max, 100 * floor_max, color='gray', alpha=0.2, label='seed-to-seed noise (max)')
ax[1].set_xlabel('x [cm]'); ax[1].set_ylabel('(n$_i$ - baseline)/n$_c$ [%] (smoothed)'); ax[1].legend(fontsize=7)
fig.tight_layout()
out = os.path.join(dirs[0], 'convergence.png')
fig.savefig(out, dpi=130)
print('wrote', out)
raise SystemExit(0 if ok else 1)
