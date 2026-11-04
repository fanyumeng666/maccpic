#!/usr/bin/env python3
"""Compare a ccp1d run with Table III of Turner et al. (2013) (data/turner2013_table3.dat).

usage: compare_table3.py <run_dir> [--case N]

Works for all four benchmark cases (Case 1 also has the full density profile,
see compare_turner.py). The case is identified from n_gas and V0 in info.txt.
Prints our value, Turner's value and the deviation for the mid-plane ion density,
the electron temperature, the electron and ion power, the ion current and the
numerical figures of merit. PASS criteria (T18): n_i within 5 %, S_e and S_i
within 10 %, J_i within 5 % of Turner's values.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import check_range, load_info, turner_table3  # noqa: E402

QE = 1.602176487e-19

args = [a for a in sys.argv[1:] if not a.startswith('--')]
d = args[0] if args else '.'
case_arg = None
for a in sys.argv[1:]:
    if a.startswith('--case'):
        case_arg = int(a.split('=', 1)[1]) if '=' in a else int(sys.argv[sys.argv.index(a) + 1])

info = load_info(os.path.join(d, 'info.txt'))
hist = np.loadtxt(os.path.join(d, 'history.dat'), ndmin=2)
here = os.path.dirname(os.path.abspath(__file__))
tab = np.loadtxt(os.path.join(here, '..', 'data', 'turner2013_table3.dat'), ndmin=2)

# identify the case
if case_arg is None:
    case, row = turner_table3(info)
    if case == 0:
        print('run parameters (n_gas = %.3e, V0 = %.0f) do not match any benchmark case' % (info['n_gas'], info['V0']))
        sys.exit(1)
else:
    case = case_arg
    row = tab[case - 1]
(_, _, _, ni_t, te_t, se_t, si_t, ji_t, wdt_t, ldx_t, nue_t, nui_t, nd_t, np_t) = row

n_avg = int(info['n_avg'])
w = hist[-n_avg:]
ne_tot, ni_tot = w[:, 2].mean(), w[:, 3].mean()
ji = 0.5 * QE * (info['Gamma_i_left'] + info['Gamma_i_right'])
lam = info['lambdaD_dx_steady'] * info['dx']
nd = info['n_e_center'] * lam / info['weight']

print('ccp1d run %s  vs  Turner et al. 2013 Table III, case %d' % (d, case))
print('%-44s %12s %12s %9s' % ('quantity', 'ccp1d', 'Turner', 'dev'))
rows = [
    ('n_i mid-plane [m^-3] (mean of 5 nodes)', info['n_i_center'], ni_t),
    ('kT_e mid-plane [eV]: 2/3 <eps_e>', info['T_eff_center'], te_t),
    ('kT_e mid-plane [eV]: m<vx^2>/e', info['Tx_e_center'], te_t),
    ('S_e electron power [W/m^2] (exact work)', info['P_part_e'], se_t),
    ('S_e electron power [W/m^2] (grid <J_e.E>)', info['P_JE_e'], se_t),
    ('S_i ion power [W/m^2] (exact work)', info['P_part_i'], si_t),
    ('J_i ion current per electrode [A/m^2]', ji, ji_t),
    ('omega_pe dt (mid-plane, steady)', info['wpe_dt_steady'], wdt_t),
    ('lambda_D / dx (mid-plane, steady)', info['lambdaD_dx_steady'], ldx_t),
    ('nu_e,max dt  (null-collision P_coll_e)', info['P_coll_e'], nue_t),
    ('nu_i,max dt  (null-collision P_coll_i)', info['P_coll_i'], nui_t),
    ('N_D particles per Debye length (mid-plane)', nd, nd_t),
    ('N_P particles: electrons', ne_tot, np_t),
    ('N_P particles: electrons + ions', ne_tot + ni_tot, np_t),
]
for name, ours, ref in rows:
    print('%-44s %12.4g %12.4g %+8.1f%%' % (name, ours, ref, 100 * (ours / ref - 1)))
print('notes: Turner does not define kT_e beyond "electron temperature in the mid-plane";')
print('       nu dt and N_P depend on the null-collision energy range and on which particles are counted;')
print('       Turner quotes 3 significant digits (n_i case 1: 1.40e14; the profile file gives 1.4046e14).')

ok = True
ok &= check_range('T18 n_i mid-plane / Turner', info['n_i_center'] / ni_t, 0.95, 1.05)
ok &= check_range('T18 S_e / Turner', info['P_part_e'] / se_t, 0.90, 1.10)
ok &= check_range('T18 S_i / Turner', info['P_part_i'] / si_t, 0.90, 1.10)
ok &= check_range('T18 J_i / Turner', ji / ji_t, 0.95, 1.05)
print('compare_table3:', 'PASS' if ok else 'FAIL')
