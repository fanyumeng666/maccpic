#!/usr/bin/env python3
"""boltzmann2term.py [out_dir] [--validate] [--bias F]: steady-state two-term
Boltzmann solver for electrons in helium with the ccp1d cross sections
(data/he), the independent reference of the swarm test T12.  Temporal-growth
(pulsed Townsend) model, as in BOLSIG+: the normalised distribution n(eps)
(per eV, integral 1) grows as exp(nu_eff t) and obeys

    nu_eff n = -d/deps [ G_E + G_el ] + S_inel
    G_E  = -D_E (n' - n/(2 eps)),  D_E = (2/3)(e/m) E^2 eps/(nu_m + nu_eff)   field heating (two-term)
    G_el = -(2m/M) nu_m eps [ n + kT (n' - n/(2 eps)) ]                       elastic cooling (Davydov)
    S_inel: excitation  -nu_k(eps) n(eps) + nu_k(eps+u_k) n(eps+u_k)
            ionization  -nu_iz(eps) n(eps) + 4 nu_iz(2eps+u) n(2eps+u)       (equal sharing, both kept)
    nu_eff = integral nu_iz n deps   (one net electron per ionization; no attachment)

with nu = N sigma v, isotropic scattering (sigma_m = sigma_elastic), zero flux at eps = 0.
Conservative Scharfetter-Gummel finite volumes for F = n/sqrt(eps); nu_eff by fixed-point iteration;
the singular system is closed by the normalisation.  Outputs: swarm_boltzmann.dat (E/N, <eps>, v_d,
nu_eff, alpha/N = nu_eff/(N v_d), power balance), swarm_b2t_eedf_<E/N>Td.dat, swarm_boltzmann_curve.dat.
--validate checks the Druyvesteyn and Maxwell limits against the analytic solutions.  --bias F
multiplies the gas density by F (the swarm's one-collision-per-step factor 0.9921)."""
import os
import sys
import numpy as np

QE, ME, MION, KB = 1.602176634e-19, 9.1093837015e-31, 6.67e-27, 1.380649e-23
here = os.path.dirname(os.path.abspath(__file__))
dhe = os.path.join(here, '..', 'data', 'he')


def load_he():
    xs = {}
    for name, thr in [('electron_scattering', False), ('excitation_1', True), ('excitation_2', True), ('ionization', True)]:
        t = np.loadtxt(os.path.join(dhe, name + '.dat'))
        xs[name] = (t[:, 0], t[:, 1], t[0, 0] if thr else 0.0, thr)
    return xs


def sigma(tab, eps):
    e, s, u, thr = tab
    out = np.interp(eps, e, s, left=0.0 if thr else s[0], right=s[-1])
    if thr:
        out = np.where(eps <= u, 0.0, out)
    return out


def bern(x):
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-8
    out[small] = 1.0 - 0.5 * x[small]
    xs = x[~small]
    out[~small] = xs / np.expm1(xs)
    return out


def solve(en_td, ngas, T_K=0.0, xs=None, h=0.1, emax=300.0, inelastic=True, growth=True, tol=1e-7, verbose=False):
    """Return dict with the grid, n(eps) and the swarm parameters (temporal-growth model if growth)."""
    if xs is None:
        xs = load_he()
    nu_eff = 0.0
    for iteration in range(30):
        r = _solve_fixed(en_td, ngas, T_K, xs, h, emax, inelastic, nu_eff if growth else 0.0, growth)
        if not growth or not inelastic:
            break
        if verbose:
            print('      iteration %d: nu_eff = %.6e' % (iteration, r['nu_iz']))
        if abs(r['nu_iz'] - nu_eff) <= tol * max(r['nu_iz'], 1.0):
            break
        nu_eff = r['nu_iz']
    r['nu_eff'] = nu_eff if growth else 0.0
    r['iterations'] = iteration + 1
    return r


def _solve_fixed(en_td, ngas, T_K, xs, h, emax, inelastic, nu_eff, growth):
    M = int(round(emax / h))
    ec = (np.arange(M) + 0.5) * h          # cell centres [eV]
    ef = np.arange(M + 1) * h              # faces
    E = en_td * 1e-21 * ngas               # V/m
    kT = KB * T_K / QE                     # eV
    vf = np.sqrt(2 * QE * ef / ME)
    sig_m = sigma(xs['electron_scattering'], ef)
    nu_m = ngas * sig_m * vf
    with np.errstate(divide='ignore', invalid='ignore'):
        D_E = (2.0 / 3.0) * (QE / ME) * E**2 * ef / (nu_m + nu_eff)      # eV^2/s
    D_E[0] = 0.0
    V_el = (2 * ME / MION) * nu_m * ef                        # eV/s
    D_T = V_el * kT
    D = D_E + D_T
    with np.errstate(divide='ignore', invalid='ignore'):
        W = np.where(ef > 0, D / (2 * ef) - V_el, 0.0)
    # unknowns n_i (per eV); fluxes written for F = n/sqrt(eps), which is smooth at eps = 0:
    #   G = -D sqrt(eps) F' - V_el sqrt(eps) F   (Scharfetter-Gummel, exact for exponential F)
    sq_c = np.sqrt(ec)
    sq_f = np.sqrt(ef)
    Dt = D * sq_f                       # diffusion of F
    Wt = -V_el * sq_f                   # drift of F (downwards)
    A = np.zeros((M, M))
    for f in range(1, M + 1):
        if Dt[f] <= 0:
            continue
        P = Wt[f] * h / Dt[f]
        aL = Dt[f] / h * bern(-P) / sq_c[f - 1]
        aR = Dt[f] / h * bern(P) / (sq_c[f] if f < M else 1.0)
        # G_f = aL n_{f-1} - aR n_f ;  cell f-1 loses G_f/h, cell f gains G_f/h ; n_M = 0
        A[f - 1, f - 1] -= aL / h
        if f < M:
            A[f - 1, f] += aR / h
            A[f, f - 1] += aL / h
            A[f, f] -= aR / h
    vc = np.sqrt(2 * QE * ec / ME)
    rates = {}
    if inelastic:
        for name in ['excitation_1', 'excitation_2', 'ionization']:
            u = xs[name][2]
            nu = ngas * sigma(xs[name], ec) * vc
            rates[name] = nu
            mult = 2.0 if (name == 'ionization' and growth) else 1.0     # electrons returned per event
            for j in np.nonzero(nu > 0)[0]:
                ed = ec[j] - u if name != 'ionization' else 0.5 * (ec[j] - u)
                if ed < 0:
                    continue
                A[j, j] -= nu[j]
                pos = ed / h - 0.5
                jd = int(np.floor(pos))
                w = pos - jd
                if jd < 0:
                    A[0, j] += mult * nu[j]
                else:
                    A[jd, j] += mult * (1 - w) * nu[j]
                    if jd + 1 < M:
                        A[jd + 1, j] += mult * w * nu[j]
    if growth and nu_eff > 0:
        A[np.arange(M), np.arange(M)] -= nu_eff          # temporal growth of the normalised distribution
    b = np.zeros(M)
    A[0, :] = h
    b[0] = 1.0
    n = np.linalg.solve(A, b)
    n = np.maximum(n, 0.0)
    n /= np.sum(n) * h
    # drift velocity from the field flux: e E v_d = integral G_E deps (per electron)
    Fc = np.concatenate([n / sq_c, [0.0]])
    GE = np.zeros(M + 1)
    for f in range(1, M + 1):
        if D_E[f] <= 0:
            continue
        GE[f] = D_E[f] * sq_f[f] / h * (Fc[f - 1] - Fc[f])
    pin = np.sum(GE) * h                      # eV/s
    vd = pin / E if E > 0 else 0.0            # m/s  (eV/s divided by V/m, eV/V = e)
    emean = np.sum(ec * n) * h
    p_el = np.sum((2 * ME / MION) * ngas * sigma(xs['electron_scattering'], ec) * vc * ec * n) * h
    p_inel = 0.0
    nu_iz = 0.0
    for name, nu in rates.items():
        u = xs[name][2]
        if name == 'ionization':
            nu_iz = np.sum(nu * n) * h
            if growth:
                p_inel += np.sum(nu * u * n) * h                      # threshold only, both electrons kept
            else:
                p_inel += np.sum(nu * 0.5 * (ec + u) * n) * h         # threshold + discarded secondary
        else:
            p_inel += np.sum(nu * u * n) * h
    p_growth = nu_eff * emean if growth else 0.0                     # energy carried by the growth (renormalisation)
    leak = GE[M] if M < len(GE) else 0.0
    return dict(ec=ec, n=n, emean=emean, vd=vd, pin=pin, p_el=p_el, p_inel=p_inel, p_growth=p_growth, nu_iz=nu_iz,
                alpha_n=nu_iz / (ngas * vd) if vd > 0 else 0.0, E=E, rates=rates, top_flux=leak)


def validate():
    ok = True
    ngas = 9.64e20
    # Druyvesteyn: constant sigma, elastic only, T = 0
    xs = load_he()
    xd = dict(xs)
    e = np.array([0.0, 1000.0]); s = np.array([5e-20, 5e-20])
    xd['electron_scattering'] = (e, s, 0.0, False)
    r = solve(10.0, ngas, 0.0, xd, h=0.05, emax=200.0, inelastic=False)
    E = 10.0 * 1e-21 * ngas
    eD = E / (ngas * 5e-20) * np.sqrt(MION / (3 * ME))       # eV
    ana = np.sqrt(r['ec']) * np.exp(-(r['ec'] / eD)**2)
    ana /= np.sum(ana) * 0.05
    m = ana > 1e-4 * ana.max()
    dev = np.max(np.abs(r['n'][m] - ana[m]) / ana.max())
    print('    Druyvesteyn (constant sigma, E/N = 10 Td, eps_D = %.2f eV): max |n - analytic| / max n = %.2e  %s' % (eD, dev, 'PASS' if dev < 1e-3 else 'FAIL'))
    ok &= dev < 1e-3
    print('        mean energy %.5f eV, analytic eps_D Gamma(5/4)/Gamma(3/4) = %.5f eV' % (r['emean'], eD * 0.906402477 / 1.225416702))
    # Maxwellian: no field, T = 1 eV, elastic only
    r = solve(1e-9, ngas, QE / KB, xd, h=0.05, emax=200.0, inelastic=False)
    print('    Maxwellian (E = 0, T = 1 eV): mean energy %.5f eV (3/2 kT = 1.5; midpoint quadrature with h = 0.05 eV)  %s' % (r['emean'], 'PASS' if abs(r['emean'] - 1.5) < 2e-3 else 'FAIL'))
    ok &= abs(r['emean'] - 1.5) < 2e-3
    # power balance of the full helium solution with growth: eEv_d = elastic + thresholds + nu_eff <eps>
    # (clean up to ~150 Td; above that the tail runs away through the top of the energy grid)
    for en, tol in ((50.0, 1e-3), (100.0, 1e-3), (150.0, 1.5e-2)):
        r = solve(en, ngas, 0.0, xs)
        bal = (r['pin'] - r['p_el'] - r['p_inel'] - r['p_growth']) / r['pin']
        print('    He %.0f Td: nu_eff = %.4e s^-1 (%d iterations), power balance (eEv_d - losses - nu_eff<eps>)/eEv_d = %.2e, top-boundary flux/nu_eff = %.1e  %s'
              % (en, r['nu_eff'], r['iterations'], bal, r['top_flux'] / max(r['nu_eff'], 1e-300), 'PASS' if abs(bal) < tol else 'FAIL'))
        ok &= abs(bal) < tol
    r = solve(200.0, ngas, 0.0, xs)
    print('    He 200 Td (information): top-boundary flux/nu_eff = %.1e, n(eps > 150 eV) = %.1e: electron runaway, the two-term'
          ' solution is no longer closed' % (r['top_flux'] / r['nu_eff'], np.sum(r['n'][r['ec'] > 150]) * 0.1))
    return ok


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    out = args[0] if args else '.'
    bias = 1.0
    if '--bias' in sys.argv:
        bias = float(sys.argv[sys.argv.index('--bias') + 1])
    if '--validate' in sys.argv:
        print('validation of the two-term solver:')
        raise SystemExit(0 if validate() else 1)
    os.makedirs(out, exist_ok=True)
    ngas = 9.64e20 * bias
    xs = load_he()
    ens = [20.0, 50.0, 100.0, 150.0, 200.0]
    with open(os.path.join(out, 'swarm_boltzmann.dat'), 'w') as f:
        f.write('# two-term Boltzmann solution (tools/boltzmann2term.py), He, n_gas = %.4e m^-3 (bias factor %.4f), T_gas = 0, temporal growth, equal energy sharing\n' % (ngas, bias))
        f.write('# 1 E/N[Td]  2 <eps>[eV]  3 v_d[m/s]  4 alpha/N = nu_eff/(N v_d)[m^2]  5 P_in=eEv_d[eV/s]  6 P_elastic  7 P_thresholds  8 nu_eff<eps>  9 balance  10 nu_eff[s^-1]\n')
        for en in ens:
            r = solve(en, ngas, 0.0, xs)
            bal = (r['pin'] - r['p_el'] - r['p_inel'] - r['p_growth']) / r['pin']
            f.write('%8.1f %12.5f %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %12.3e %14.6e\n' % (en, r['emean'], r['vd'], r['alpha_n'], r['pin'], r['p_el'], r['p_inel'], r['p_growth'], bal, r['nu_eff']))
            with open(os.path.join(out, 'swarm_b2t_eedf_%dTd.dat' % int(en)), 'w') as g:
                g.write('# two-term EEDF n(eps) at E/N = %.0f Td, integral 1\n# 1 eps[eV]  2 n[eV^-1]\n' % en)
                for e_, n_ in zip(r['ec'], r['n']):
                    g.write('%12.4f %15.7e\n' % (e_, n_))
            print('E/N = %5.0f Td: <eps> = %7.3f eV  v_d = %.4e m/s  nu_eff = %.3e s^-1  alpha/N = %.3e m^2  balance %.1e  top flux/nu_eff %.1e' % (en, r['emean'], r['vd'], r['nu_eff'], r['alpha_n'], bal, r['top_flux'] / max(r['nu_eff'], 1e-300)))
    with open(os.path.join(out, 'swarm_boltzmann_curve.dat'), 'w') as f:
        f.write('# two-term Boltzmann: 1 E/N[Td] 2 <eps>[eV] 3 v_d[m/s] 4 alpha/N[m^2]\n')
        for en in np.logspace(np.log10(10), np.log10(300), 13):
            r = solve(en, ngas, 0.0, xs)
            f.write('%10.3f %12.5f %14.6e %14.6e\n' % (en, r['emean'], r['vd'], r['alpha_n']))
    print('wrote', os.path.join(out, 'swarm_boltzmann.dat'))
