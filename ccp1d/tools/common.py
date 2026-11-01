"""Shared helpers for the ccp1d analysis scripts (numpy + matplotlib only)."""
import os
import re
import sys
import numpy as np

QE = 1.602176634e-19
ME = 9.1093837015e-31
MION = 6.67e-27
EPS0 = 8.8541878128e-12
KB = 1.380649e-23


def load_dat(path):
    """Return (array, header_lines) of a '#'-commented ASCII table."""
    header = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                header.append(line[1:].strip())
            else:
                break
    return np.loadtxt(path, comments='#'), header


def load_info(path):
    """Parse info.txt ('key = value  # unit') into a dict; numbers become floats."""
    info = {}
    with open(path) as f:
        for line in f:
            line = line.split('#', 1)[0].strip()
            if '=' not in line:
                continue
            key, val = [s.strip() for s in line.split('=', 1)]
            try:
                info[key] = float(val.replace('d', 'e').replace('D', 'e'))
            except ValueError:
                info[key] = val
    return info


def run_dir_from_argv(default='.'):
    d = sys.argv[1] if len(sys.argv) > 1 else default
    if not os.path.isdir(d):
        sys.exit('run directory not found: ' + d)
    return d


def history_columns():
    """Column index (0-based) by name for history.dat."""
    names = ['period', 't', 'N_e', 'N_i', 'n_e_mean', 'n_i_mean', 'n_e_center', 'n_i_center',
             'Gi_L', 'Gi_R', 'Ge_L', 'Ge_R', 'Ei_L', 'Ei_R', 'P_e', 'P_i', 'P_VI',
             'n_ioniz', 'n_exc', 'KE_e', 'KE_i', 'W_field', 'E_wall_cum', 'Work_cum',
             'gauss_res_max', 'wall', 'Tx_e', 'Tperp_e', 'Tx_i', 'Tperp_i', 'numax_exceed']
    return {n: i for i, n in enumerate(names)}


def profile_columns():
    names = ['x', 'n_e', 'n_i', 'phi', 'E', 'eps_e', 'ioniz', 'PeE', 'PiE', 'flux_i', 'eps_x_e', 'eps_i']
    return {n: i for i, n in enumerate(names)}


def check(name, value, tol, fmt='%.3e'):
    ok = value <= tol
    print('    %s %s = %s  (tol %s)' % ('PASS' if ok else 'FAIL', name, fmt % value, fmt % tol))
    return ok


def check_range(name, value, lo, hi, fmt='%.3e'):
    ok = lo <= value <= hi
    print('    %s %s = %s  (range %s .. %s)' % ('PASS' if ok else 'FAIL', name, fmt % value, fmt % lo, fmt % hi))
    return ok
