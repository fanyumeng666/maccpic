#!/usr/bin/env python3
"""compare_swarm.py <dir>: T12.  Compares the PIC swarm (swarm_pic.dat, swarm_eedf_*.dat)
with the two-term Boltzmann solution (swarm_boltzmann.dat, swarm_b2t_eedf_*.dat) at
E/N = 50, 100, 200, 500 Td: mean energy, drift velocity, alpha/N, energy balance.
Tolerances: 3 % for <eps> and v_d up to 200 Td, 8 % at 500 Td (two-term approximation),
15 % for alpha/N (tail sensitive); the PIC energy balance must close.  Writes swarm.png."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import load_dat, check

d = sys.argv[1] if len(sys.argv) > 1 else '.'
pic, hp = load_dat(os.path.join(d, 'swarm_pic.dat'))
b2t, hb = load_dat(os.path.join(d, 'swarm_boltzmann.dat'))
curve, _ = load_dat(os.path.join(d, 'swarm_boltzmann_curve.dat'))
ok = True
print('%6s | %9s %9s %7s | %10s %10s %7s | %10s %10s %7s | %9s %9s' % ('E/N', '<eps>PIC', '<eps>B2T', 'dev', 'v_d PIC', 'v_d B2T', 'dev', 'a/N PIC', 'a/N B2T', 'dev', 'bal PIC', 'ident'))
for i in range(len(pic)):
    en = pic[i, 0]
    j = np.argmin(np.abs(b2t[:, 0] - en))
    de = pic[i, 2] / b2t[j, 1] - 1
    dv = pic[i, 3] / b2t[j, 2] - 1
    da = pic[i, 4] / b2t[j, 3] - 1 if b2t[j, 3] > 0 else 0.0
    print('%6.0f | %9.3f %9.3f %+6.1f%% | %10.3e %10.3e %+6.1f%% | %10.3e %10.3e %+6.1f%% | %9.1e %9.1e' %
          (en, pic[i, 2], b2t[j, 1], 100 * de, pic[i, 3], b2t[j, 2], 100 * dv, pic[i, 4], b2t[j, 3], 100 * da, pic[i, 7], pic[i, 12]))
    ok &= check('T12 %3.0f Td PIC exact identity work - dKE - losses (relative)' % en, abs(pic[i, 12]), 1e-9)
    if en > 150:
        print('    INFO %3.0f Td: electron runaway regime (friction force < eE above ~100 eV); PIC and two-term are compared for information only' % en)
        continue
    tol = 0.03 if en <= 100 else 0.05
    ok &= check('T12 %3.0f Td |<eps> PIC/B2T - 1|' % en, abs(de), tol)
    ok &= check('T12 %3.0f Td |v_d PIC/B2T - 1|' % en, abs(dv), tol)
    if b2t[j, 3] > 0 and pic[i, 4] > 0:
        ok &= check('T12 %3.0f Td |alpha/N PIC/B2T - 1|' % en, abs(da), 0.15)
    ok &= check('T12 %3.0f Td PIC stationary balance |eEv_d - losses|/eEv_d' % en, abs(pic[i, 7]), 0.01)
print('T12:', 'PASS' if ok else 'FAIL')

fig, ax = plt.subplots(1, 3, figsize=(17, 5))
a = ax[0]
a.loglog(curve[:, 0], curve[:, 1], 'k-', label='two-term Boltzmann')
a.loglog(pic[:, 0], pic[:, 2], 'o', color='C3', label='PIC swarm (ccp1d push + mcc)')
a.set_xlabel('E/N [Td]'); a.set_ylabel('mean energy [eV]'); a.legend(fontsize=8); a.set_title('mean electron energy', fontsize=10)
a = ax[1]
a.loglog(curve[:, 0], curve[:, 2], 'k-', label='two-term Boltzmann')
a.loglog(pic[:, 0], pic[:, 3], 'o', color='C3', label='PIC swarm')
a.set_xlabel('E/N [Td]'); a.set_ylabel('drift velocity [m/s]'); a.legend(fontsize=8); a.set_title('drift velocity <v$_x$>', fontsize=10)
a2 = a.twinx()
m = curve[:, 3] > 0
a2.loglog(curve[m, 0], curve[m, 3], 'k--', lw=0.8)
mp = pic[:, 4] > 0
a2.loglog(pic[mp, 0], pic[mp, 4], 's', color='C0', label='alpha/N (right axis)')
a2.set_ylabel('alpha/N [m$^2$]'); a2.legend(fontsize=8, loc='lower right')
a = ax[2]
for i, en in enumerate(pic[:, 0]):
    pe, _ = load_dat(os.path.join(d, 'swarm_eedf_%dTd.dat' % int(en)))
    be, _ = load_dat(os.path.join(d, 'swarm_b2t_eedf_%dTd.dat' % int(en)))
    m = pe[:, 1] > 0
    a.semilogy(pe[m, 0], pe[m, 1], '.', ms=2.5, color='C%d' % i, label='PIC %.0f Td' % en)
    a.semilogy(be[:, 0], be[:, 1], '-', lw=0.9, color='C%d' % i)
a.set_ylim(1e-7, 1); a.set_xlim(0, 200); a.set_xlabel('electron energy [eV]'); a.set_ylabel('n($\\epsilon$) [eV$^{-1}$]')
a.set_title('energy distributions: PIC (dots) vs two-term (lines)', fontsize=10); a.legend(fontsize=8)
fig.suptitle('T12 electron swarm in He at n$_g$ = 9.64e20 m$^{-3}$ (pulsed Townsend: secondaries kept, population renormalised)', fontsize=10)
fig.tight_layout()
out = os.path.join(d, 'swarm.png')
fig.savefig(out, dpi=130)
print('wrote', out)
raise SystemExit(0 if ok else 1)
