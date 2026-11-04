#!/usr/bin/env python3
"""compare_cases.py <run_dir> [<run_dir> ...] [--out file.png]

Pressure dependence across Turner benchmark cases, drawn in the units of the figures of
Turner et al. (2013): ion density (Figs 7-10, 1e15 m^-3), ionization source (Fig. 3,
1e20 m^-3 s^-1), power density coupled to electrons and ions (Figs 4-5, kW/m^3), the
EEPF normalised to int sqrt(eps) f deps = 1 (Fig. 6, eV^-3/2) and the IEDF at the driven
electrode. Stars mark the mid-plane densities of Turner's Table III. Also prints the
Table III quantities of every run. Default output: results/case_comparison.png."""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_dat, load_info, profile_columns, turner_table3, QE  # noqa: E402

KB = 1.3806504e-23
argv = sys.argv[1:]
out = None
if '--out' in argv:
    i = argv.index('--out'); out = argv[i + 1]; del argv[i:i + 2]
dirs = argv
if not dirs:
    sys.exit(__doc__)
if out is None:
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'case_comparison.png')

pc = profile_columns()
fig, ax = plt.subplots(2, 3, figsize=(16, 9))
print('%-28s %6s %5s %10s %7s %7s %7s %7s' % ('run', 'mTorr', 'V0', 'n_i,c', 'T_eff', 'S_e', 'S_i', 'J_i'))
emax = 0.0
for k, d in enumerate(dirs):
    info = load_info(os.path.join(d, 'info.txt'))
    prof, _ = load_dat(os.path.join(d, 'profiles.dat'))
    eepf, _ = load_dat(os.path.join(d, 'eepf.dat'))
    iedf, _ = load_dat(os.path.join(d, 'iedf_left.dat'))
    case, row = turner_table3(info)
    p_mtorr = info['n_gas'] * KB * info['T_gas'] / 133.322 * 1e3
    lab = ('case %d: ' % case if case else '') + '%.0f mTorr, %.0f V' % (p_mtorr, info['V0'])
    c = 'C%d' % k
    xl = prof[:, pc['x']] / info['L']
    ji = 0.5 * QE * (info['Gamma_i_left'] + info['Gamma_i_right'])
    print('%-28s %6.0f %5.0f %10.3e %7.2f %7.1f %7.1f %7.3f' % (d, p_mtorr, info['V0'], info['n_i_center'], info['T_eff_center'], info['P_part_e'], info['P_part_i'], ji))
    if row is not None:
        print('%-28s %6s %5s %10.3e %7.2f %7.1f %7.1f %7.3f' % ('   Turner Table III', '', '', row[3], row[4], row[5], row[6], row[7]))
    a = ax[0, 0]; a.plot(xl, prof[:, pc['n_i']] * 1e-15, color=c, label=lab)
    if row is not None:
        a.plot(0.5, row[3] * 1e-15, '*', color=c, ms=12, mec='k', label='Turner Table III mid-plane' if k == 0 else None)
    a = ax[0, 1]; a.plot(xl, prof[:, pc['ioniz']] * 1e-20, color=c, label=lab)
    a = ax[0, 2]; a.plot(xl, prof[:, pc['P_work_e']] * 1e-3, color=c, label=lab + ' (exact work)')
    a.plot(xl, prof[:, pc['PeE']] * 1e-3, ':', color=c, lw=1, label='grid <J$_e$.E>' if k == 0 else None)
    a = ax[1, 0]; a.plot(xl, prof[:, pc['P_work_i']] * 1e-3, color=c, label=lab)
    e, F, fp = eepf[:, 0], eepf[:, 1], eepf[:, 2]
    ne_c = np.sum(F) * (e[1] - e[0])
    m = fp > 0
    a = ax[1, 1]; a.semilogy(e[m], fp[m] / ne_c, '.', ms=3, color=c, label=lab + ' (T$_{eff}$ = %.2f eV)' % info['T_eff_center'])
    ei, fi = iedf[:, 0], iedf[:, 1]
    tot = np.sum(fi) * (ei[1] - ei[0])
    emax = max(emax, ei[np.nonzero(fi)[0][-1]] if np.any(fi > 0) else 0.0)
    a = ax[1, 2]; a.plot(ei, fi / tot, color=c, label=lab + ' (mean %.0f eV)' % info['E_i_mean_left'])
ax[0, 0].set_ylabel('n$_i$ [10$^{15}$ m$^{-3}$]'); ax[0, 0].set_title('time-averaged ion density (Turner Figs 7-10)'); ax[0, 0].set_ylim(bottom=0)
ax[0, 1].set_ylabel('R$_i$ [10$^{20}$ m$^{-3}$ s$^{-1}$]'); ax[0, 1].set_title('ionization source (Turner Fig. 3)'); ax[0, 1].set_ylim(bottom=0)
ax[0, 2].set_ylabel('P$_e$ [kW/m$^3$]'); ax[0, 2].set_title('power density coupled to electrons (Turner Fig. 4)'); ax[0, 2].axhline(0, color='k', lw=0.6)
ax[1, 0].set_ylabel('P$_i$ [kW/m$^3$]'); ax[1, 0].set_title('power density coupled to ions (Turner Fig. 5)'); ax[1, 0].set_ylim(bottom=0)
ax[1, 1].set_xlabel(r'$\epsilon$ [eV]'); ax[1, 1].set_ylabel(r'f($\epsilon$) [eV$^{-3/2}$]'); ax[1, 1].set_xlim(0, 60); ax[1, 1].set_ylim(1e-5, 1e-1)
ax[1, 1].set_title(r'EEPF, central 10 %, $\int\sqrt{\epsilon}\,f\,d\epsilon$ = 1 (Turner Fig. 6)')
ax[1, 2].set_xlabel(r'$\epsilon_i$ [eV]'); ax[1, 2].set_ylabel('normalised IEDF [eV$^{-1}$]'); ax[1, 2].set_title('ion energy distribution at the driven electrode'); ax[1, 2].set_ylim(bottom=0); ax[1, 2].set_xlim(0, 1.1 * emax if emax > 0 else None)
for a in (ax[0, 0], ax[0, 1], ax[0, 2], ax[1, 0]):
    a.set_xlabel('x / L'); a.set_xlim(0, 1)
for row_ in ax:
    for a in row_:
        a.legend(fontsize=8)
fig.suptitle('ccp1d: pressure dependence across the Turner et al. (2013) benchmark cases', fontsize=12)
fig.tight_layout()
fig.savefig(out, dpi=120)
print('wrote', out)
