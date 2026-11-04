#!/usr/bin/env python3
"""theory_check.py <run_dir>: compare a ccp1d run with classical low-temperature plasma theory.

Every line prints the simulation value, the theoretical estimate, their ratio and the formula.
Sources: Lieberman & Lichtenberg, Principles of Plasma Discharges (2nd ed., "L&L"), Godyak 1986,
Turner 2013. Needs info.txt, profiles.dat, eepf.dat, iedf_left.dat and (for the sheath motion
and the current amplitude) xt_ne.dat, xt_ni.dat, xt_Je.dat of the run directory."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_dat, load_info, profile_columns, QE  # noqa: E402

ME, MI, EPS0, KB = 9.10938215e-31, 6.67e-27, 8.854187817e-12, 1.3806504e-23
d = sys.argv[1] if len(sys.argv) > 1 else '.'
info = load_info(os.path.join(d, 'info.txt'))
prof, _ = load_dat(os.path.join(d, 'profiles.dat'))
eepf, _ = load_dat(os.path.join(d, 'eepf.dat'))
iedf, _ = load_dat(os.path.join(d, 'iedf_left.dat'))
pc = profile_columns()
here = os.path.dirname(os.path.abspath(__file__))
xs = {k: np.loadtxt(os.path.join(here, '..', 'data', 'he', k + '.dat')) for k in
      ['electron_scattering', 'excitation_1', 'excitation_2', 'ionization', 'ion_scattering', 'ion_back_scatter']}


def sig(name, e):
    t = xs[name]
    return np.interp(e, t[:, 0], t[:, 1], left=0.0 if name != 'ion_scattering' and name != 'ion_back_scatter' else t[0, 1], right=t[-1, 1])


L, ng, V0, f = info['L'], info['n_gas'], info['V0'], info['freq']
w = 2 * np.pi * f
x = prof[:, pc['x']]; ne = prof[:, pc['n_e']]; ni = prof[:, pc['n_i']]; phi = prof[:, pc['phi']]
flux = prof[:, pc['flux_i']]; rate = prof[:, pc['ioniz']]; pwe = prof[:, pc['P_work_e']]
nx = len(x) - 1; dx = x[1] - x[0]; mid = nx // 2
Te = info['T_eff_center']; Tx = info['Tx_e_center']; n0 = info['n_i_center']; ne0 = info['n_e_center']
Gi = info['Gamma_i_left'] + info['Gamma_i_right']; Ge = info['Gamma_e_left'] + info['Gamma_e_right']
uB = np.sqrt(QE * Te / MI)
vbar = np.sqrt(8 * QE * Te / (np.pi * ME))
phic = phi.max()

# EEPF of the central 10 %: F [m^-3 eV^-1]; rate coefficients <sigma v>
e_b, F = eepf[:, 0], eepf[:, 1]
de = e_b[1] - e_b[0]
nc = F.sum() * de
v_b = np.sqrt(2 * QE * e_b / ME)


def K_eepf(name):
    return np.sum(sig(name, e_b) * v_b * F) * de / nc


def K_maxw(name, T, loss=False):
    e = np.linspace(0.005, 400, 80000)
    FM = 2 * np.sqrt(e / np.pi) * T ** -1.5 * np.exp(-e / T)
    v = np.sqrt(2 * QE * e / ME)
    g = sig(name, e) * v * FM
    if loss:
        g = g * 2 * ME / MI * e
    return np.trapezoid(g, e)


sig_tot_e = lambda e: sig('electron_scattering', e) + sig('excitation_1', e) + sig('excitation_2', e) + sig('ionization', e)
nu_m = ng * np.sum(sig_tot_e(e_b) * v_b * F) * de / nc
eps_mean = np.sum(e_b * F) * de / nc
lam_m = 1.0 / (ng * sig_tot_e(eps_mean))
lam_inel = 1.0 / (ng * (sig('excitation_1', 30.0) + sig('excitation_2', 30.0) + sig('ionization', 30.0)))
lam_eps = lam_m * np.sqrt(MI / (2 * ME))
T_lo = T_hi = None
for lo, hi in ((1, 8), (25, 80)):
    m_ = (e_b >= lo) & (e_b <= hi) & (F > 0)
    if m_.sum() > 3:
        T_ = -1.0 / np.polyfit(e_b[m_], np.log(F[m_] / np.sqrt(e_b[m_])), 1)[0]
        if lo == 1:
            T_lo = T_
        else:
            T_hi = T_

# sheath motion and current amplitude from the x-t maps (if present)
have_xt = all(os.path.exists(os.path.join(d, 'xt_%s.dat' % k)) for k in ('ne', 'ni', 'Je'))
if have_xt:
    xne, _ = load_dat(os.path.join(d, 'xt_ne.dat')); xni, _ = load_dat(os.path.join(d, 'xt_ni.dat')); xje, _ = load_dat(os.path.join(d, 'xt_Je.dat'))
    sl = []; sr = []
    for p in range(xne.shape[0]):
        r = xne[p] / np.maximum(xni[p], 1e-300)
        idx = np.where(r >= 0.5)[0]
        sl.append(x[idx[0]]); sr.append(L - x[idx[-1]])
    sl = np.array(sl); sr = np.array(sr)
    s_max = 0.5 * (sl.max() + sr.max()); s_min = 0.5 * (sl.min() + sr.min())
    nph = xje.shape[0]
    jc = xje[:, mid - 2:mid + 3].mean(axis=1)
    J0 = 2.0 / nph * abs(np.sum(jc * np.exp(-2j * np.pi * np.arange(nph) / nph)))
else:
    s_max = s_min = J0 = None

print('theory_check: %s   (Turner case: %s mTorr, %.0f V; T_eff = %.2f eV, T_x = %.2f eV, n_i,c = %.3e)' %
      (d, '%.0f' % (ng * KB * info['T_gas'] / 133.322e-3), V0, Te, Tx, n0))
print('regime: nu_m/omega = %.2f (nu_m = %.2e s^-1 from the central EEPF, sigma_m = all isotropic processes)' % (nu_m / w, nu_m))
print('        electron mean free path at <eps> = %.1f eV: %.2f cm = %.2f L; energy relaxation length (elastic) %.0f L;' % (eps_mean, lam_m * 1e2, lam_m / L, lam_eps / L))
print('        inelastic mean free path at 30 eV: %.2f L; ion mean free path at T_e/4: %.2f mm' % (lam_inel / L, 1e3 / (ng * (sig('ion_scattering', Te / 4) + sig('ion_back_scatter', Te / 4)))))
if T_lo and T_hi:
    print('        EEPF slope temperatures: 1-8 eV %.2f eV, 25-80 eV %.2f eV' % (T_lo, T_hi))
print()
rows = []


def row(name, sim, th, formula):
    rows.append((name, sim, th, formula))


# 1 plasma potential
row('time-averaged plasma potential [V]', phic, V0 / 2, 'symmetric discharge: <phi_c> = V0/2 (each sheath swings 0..V0)')
# 2 ion energies
emax = iedf[np.nonzero(iedf[:, 1])[0][-1], 0]
row('maximum ion energy at the electrode [eV]', emax, phic, 'ions see the time-averaged sheath drop e<phi_c> when omega tau_i >> 1')
row('field work per ion S_i/(e Gamma_i) [eV]', info['P_part_i'] / (QE * Gi), np.sum(rate * phi) / np.sum(rate),
    'ions born at rest at the ionization-weighted potential <phi_birth>, all work ends at the wall or in the gas')
row('ion energy to the wall per ion [eV]', info['P_wall_i'] / (QE * Gi), info['E_i_mean_left'], 'consistency with the IEDF mean (energy lost to gas by charge exchange is the rest)')
row('electron energy to the wall per electron [eV]', info['P_wall_e'] / (QE * Ge), 2 * Te, 'Maxwellian flux through a retarding sheath carries 2 T_e per electron (L&L 2.4)')
if T_hi:
    row('   same, with the hot-group temperature', info['P_wall_e'] / (QE * Ge), 2 * T_hi, '2 T_hot with T_hot = %.2f eV from the EEPF slope 25-80 eV: only the tail reaches the wall' % T_hi)
# 3 sheath: Child law
if have_xt:
    Ji = 0.5 * QE * Gi
    s_cl = np.sqrt(200.0 / 243.0 * EPS0 * np.sqrt(2 * QE / MI) * phic ** 1.5 / Ji)
    lam_i_s = 1.0 / (ng * (sig('ion_scattering', phic / 4) + sig('ion_back_scatter', phic / 4)))
    s_col = (1.68 * EPS0 * np.sqrt(2 * QE * lam_i_s / (np.pi * MI)) * phic ** 1.5 / Ji) ** 0.4
    row('maximum sheath thickness [cm]', s_max * 1e2, s_cl * 1e2, 'collisionless RF Child law s_m = [200/243 eps0 (2e/M)^1/2 Vbar^3/2 / J_i]^1/2, Vbar = <phi_c> (L&L 11.2)')
    row('   same, collisional (DC) Child law [cm]', s_max * 1e2, s_col * 1e2, 's = [1.68 eps0 (2 e lambda_i / pi M)^1/2 V^3/2 / J_i]^2/5, lambda_i = %.1f mm at eps_rel = Vbar/4 (L&L 6.5); s_m/lambda_i = %.0f' % (lam_i_s * 1e3, s_max / lam_i_s))
    tau_i = 3 * s_max / np.sqrt(2 * QE * phic / MI)
    row('ion sheath transit time / RF period', tau_i * f, 1.0, 'tau_i = 3 s_m / sqrt(2 e Vbar / M) (Child-law sheath); >> 1 means the ions average the RF field')
# 4 Bohm point and edge-to-centre density ratio
vi = np.abs(flux) / np.maximum(ni, 1e-300)


def bohm_point(T):
    ib = np.where(vi[:mid] <= np.sqrt(QE * T / MI))[0]
    return x[ib[0]] if len(ib) else np.nan


if have_xt:
    i_avg = np.where(ne[:mid] / np.maximum(ni[:mid], 1e-300) >= 0.5)[0][0]
    row('Bohm point (|v_i| = u_B(T_eff)) [cm]', bohm_point(Te) * 1e2, s_max * 1e2, 'ions must reach u_B = (e T_e/M)^1/2 before entering the sheath (Bohm criterion, L&L 6.2); compared with the maximum sheath extension; the time-averaged edge n_e/n_i = 0.5 is at %.2f cm' % (x[i_avg] * 1e2))
    if T_lo:
        row('   same with u_B of the cold group [cm]', bohm_point(T_lo) * 1e2, s_max * 1e2, 'u_B with T_e = %.2f eV (the low-energy electrons set the presheath, L&L 6.2 / two-temperature Bohm speed)' % T_lo)
h_meas = 0.5 * Gi / (n0 * uB)
lam_i = 1.0 / (ng * (sig('ion_scattering', Te / 4) + sig('ion_back_scatter', Te / 4)))
h_l = 0.86 / np.sqrt(3 + L / (2 * lam_i))
row('edge/centre density factor h_l = Gamma_i/(n_0 u_B)', h_meas, h_l, 'Godyak: h_l = 0.86 (3 + L/2 lambda_i)^-1/2, lambda_i = %.1f mm (L&L 10.2); assumes a thin sheath' % (lam_i * 1e3))
if have_xt:
    lb = L - 2 * s_max
    row('   same with the never-sheath bulk length', h_meas, 0.86 / np.sqrt(3 + lb / (2 * lam_i)), 'plasma length l = L - 2 s_m = %.2f cm between the maximum sheath positions' % (lb * 1e2))
# 5 Boltzmann relation in the bulk
sel = (phic - phi < 2 * Te) & (ne > 0)
if have_xt:
    sel &= (x > s_max) & (x < L - s_max)
if sel.sum() > 5:
    TB = -1.0 / np.polyfit(phic - phi[sel], np.log(ne[sel]), 1)[0]
    row('Boltzmann temperature of the bulk n_e(phi) [eV]', TB, Te, 'n_e = n_0 exp(-(phi_0 - phi)/T_B) fitted over %d nodes never covered by a sheath with phi_0 - phi < 2 T_eff (drop %.1f V); theory T_B = T_e%s' %
        (sel.sum(), (phic - phi[sel]).max(), ' (cold group %.2f eV)' % T_lo if T_lo else ''))
# 6 ohmic heating at the centre
if have_xt:
    P_ohm = J0 ** 2 * ME * nu_m / (2 * ne0 * QE ** 2)
    P_sim = pwe[int(0.4 * nx):int(0.6 * nx) + 1].mean()
    row('electron heating at the centre [kW/m^3]', P_sim * 1e-3, P_ohm * 1e-3, 'Drude ohmic heating J0^2 m nu_m / (2 n e^2) with J0 = %.1f A/m^2 (fundamental of J_e at the centre); local theory, needs nu_m >~ omega and lambda_m << L' % J0)
# 7 hard-wall stochastic heating
if have_xt:
    u0 = w * 0.5 * (s_max - s_min)
    i_s = int(round(0.5 * (s_max + s_min) / dx))
    n_s = ni[i_s]
    S_stoc = 0.5 * ME * vbar * n_s * u0 ** 2
    row('electron power S_e [W/m^2]', info['P_part_e'], 2 * S_stoc, 'hard-wall (Fermi) stochastic heating 2 x 1/2 m vbar_e n_s u0^2, u0 = omega x half the sheath excursion = %.2e m/s, n_s = %.2e (L&L 11.2); upper estimate' % (u0, n_s))
    bulk = (x > s_max) & (x < L - s_max) & (ne > 0)
    S_ohm = np.trapezoid(J0 ** 2 * ME * nu_m / (2 * ne[bulk] * QE ** 2), x[bulk])
    row('   same, stochastic + bulk ohmic', info['P_part_e'], 2 * S_stoc + S_ohm, 'adds the Drude ohmic heating integrated over the never-sheath bulk (%.1f W/m^2)' % S_ohm)
# 8 energy cost per electron-ion pair
eps_c_sim = (info['P_inelastic'] + info['P_elastic_e']) / (QE * Gi)
Km = {k: K_maxw(k, Te) for k in ['excitation_1', 'excitation_2', 'ionization']}
eps_c_M = (19.82 * Km['excitation_1'] + 20.61 * Km['excitation_2'] + 24.59 * Km['ionization'] + K_maxw('electron_scattering', Te, loss=True)) / Km['ionization']
Ke = {k: K_eepf(k) for k in ['excitation_1', 'excitation_2', 'ionization']}
Kel_loss = np.sum(sig('electron_scattering', e_b) * v_b * F * 2 * ME / MI * e_b) * de / nc
eps_c_E = (19.82 * Ke['excitation_1'] + 20.61 * Ke['excitation_2'] + 24.59 * Ke['ionization'] + Kel_loss) / Ke['ionization']
row('collisional energy cost per e-i pair eps_c [eV]', eps_c_sim, eps_c_M, 'eps_c = sum_k eps_k K_k / K_iz with Maxwellian rates at T_eff (L&L 3.5); simulation: (P_inel + P_el)/(e Gamma_i)')
row('   same with the central EEPF instead of a Maxwellian [eV]', eps_c_sim, eps_c_E, 'rates integrated over the measured central EEPF (checks the non-Maxwellian tail)')
# 9 particle balance
Ne_line = np.trapezoid(ne, x)
R_sim = Gi
row('ionization rate per area [m^-2 s^-1]', R_sim, ng * Km['ionization'] * Ne_line, 'particle balance n_g K_iz(T_eff, Maxwellian) int n_e dx = Gamma_i,left + Gamma_i,right (global model, L&L 10.2)')
row('   same with the central EEPF', R_sim, ng * Ke['ionization'] * Ne_line, 'n_g K_iz(EEPF) int n_e dx; ratio < 1 means the tail is weaker away from the centre, > 1 stronger')
lo_, hi_ = 1.0, 40.0
for _ in range(60):
    T_gm = 0.5 * (lo_ + hi_)
    if ng * K_maxw('ionization', T_gm) * Ne_line < Gi:
        lo_ = T_gm
    else:
        hi_ = T_gm
row('electron temperature T_eff [eV]', Te, T_gm, 'T_e that balances ionization and ion loss with Maxwellian rates (global model); the real EEPF is not Maxwellian')
# 10 sheath capacitance / current amplitude
if have_xt:
    row('RF current amplitude J0 [A/m^2]', J0, 1.52 * EPS0 * w * V0 / (2 * s_max), 'two Child-law sheaths in series, C_s = 1.52 eps0 A/s_m (L&L 11.2): J0 = omega C_s V0/2 A; order of magnitude')

print('%-52s %11s %11s %7s' % ('quantity', 'simulation', 'theory', 'ratio'))
for name, s, t, formula in rows:
    print('%-52s %11.4g %11.4g %7.3f   %s' % (name, s, t, s / t if t else np.nan, formula))
