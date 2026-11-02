#!/usr/bin/env python3
"""plot_power.py <run_dir>: where the absorbed power goes (from the energy
balance in info.txt) -> power.png and a printed table."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import load_info, run_dir_from_argv

d = run_dir_from_argv()
i = load_info(os.path.join(d, 'info.txt'))
src = [('source <V I>', i['P_VI']), ('work on electrons', i['P_part_e']), ('work on ions', i['P_part_i'])]
sinks = [('ions to electrodes', i['P_wall_i']), ('electrons to electrodes', i['P_wall_e']), ('excitation + ionization', i['P_inelastic']),
         ('ions to gas (charge exchange, elastic)', i['P_gas_i']), ('electron elastic loss', i['P_elastic_e']),
         ('stored energy change', i['P_dKE_e'] + i['P_dKE_i'] + i['dW_field_dt'])]
tot = i['P_part_total'] + i['P_newion']
print('power budget [W/m^2] (window average):')
for n, v in src:
    print('  %-42s %8.2f' % (n, v))
print('  sinks:')
for n, v in sinks:
    print('  %-42s %8.2f  (%5.1f %%)' % (n, v, 100 * v / tot))
print('  %-42s %8.2e' % ('closure residual', i['closure_residual']))
print('  %-42s %8.2f  (%.2f %% of P_part)' % ('numerical energy error P_VI - dW/dt - P_part', i['P_numerical'], 100 * i['P_numerical_relative']))
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
ax[0].barh([n for n, _ in src], [v for _, v in src], color=['k', 'C0', 'C1'])
ax[0].set_xlabel('W/m$^2$'); ax[0].set_title('power in: source and work on each species'); ax[0].invert_yaxis()
for k, (n, v) in enumerate(src):
    ax[0].text(v, k, ' %.1f' % v, va='center')
ax[1].barh([n for n, _ in sinks], [v for _, v in sinks], color='C3')
ax[1].set_xlabel('W/m$^2$'); ax[1].set_title('power out (sum = %.1f, residual %.0e)' % (sum(v for _, v in sinks), i['closure_residual'])); ax[1].invert_yaxis()
for k, (n, v) in enumerate(sinks):
    ax[1].text(max(v, 0), k, ' %.1f (%.0f %%)' % (v, 100 * v / tot), va='center')
fig.suptitle(os.path.abspath(d), fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(d, 'power.png'), dpi=120)
print('wrote', os.path.join(d, 'power.png'))
