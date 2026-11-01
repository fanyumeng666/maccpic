#!/usr/bin/env python3
"""plot_initial.py <run_dir> [out_dir]: figures of the initial state of a run
made with dump_initial = .true. (init_particles_rank*.dat, init_grid.dat,
info.txt) plus the collision set-up from data/he:
  init_phase_space.png   x-vx phase space of electrons and ions at t = 0
  init_velocities.png    vx histograms and the electron energy distribution vs Maxwellian
  init_grid.png          n_e, n_i, rho, phi, E at t = 0 (shot noise of the loading)
  init_drive.png         V(t) over one period with the time steps, and the scale comparison
  init_xsec.png          cross sections, collision frequencies, nu_max and the null-collision probability
"""
import glob
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import load_dat, load_info, run_dir_from_argv, QE, ME, MION, EPS0, KB

d = run_dir_from_argv()
out = sys.argv[2] if len(sys.argv) > 2 else d
os.makedirs(out, exist_ok=True)
info = load_info(os.path.join(d, 'info.txt'))
L, nx, dx, dt = info['L'], int(info['nx']), info['dx'], info['dt']
n0, Te0, Ti0, w = info['n0'], info['Te0'], info['Ti0'], info['weight']
V0, freq, ngas = info['V0'], info['freq'], info['n_gas']
files = sorted(glob.glob(os.path.join(d, 'init_particles_rank*.dat')))
parts = np.vstack([np.loadtxt(f, comments='#') for f in files])
el = parts[parts[:, 0] == 1]
io = parts[parts[:, 0] == 2]
grid, _ = load_dat(os.path.join(d, 'init_grid.dat'))
xg = grid[:, 0]
print('particles: %d electrons, %d ions on %d rank files' % (len(el), len(io), len(files)))

vte = np.sqrt(KB * Te0 / ME)
vti = np.sqrt(KB * Ti0 / MION)
te_ev = KB * Te0 / QE
lam_d = np.sqrt(EPS0 * KB * Te0 / (n0 * QE**2))
wpe = np.sqrt(n0 * QE**2 / (EPS0 * ME))

# ---- 1. phase space ----
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
sub = np.random.default_rng(0).choice(len(el), size=min(20000, len(el)), replace=False)
ax[0].plot(el[sub, 1] * 1e2, el[sub, 2] * 1e-6, ',', color='C0', alpha=0.6)
ax[0].set_xlabel('x [cm]'); ax[0].set_ylabel('v$_x$ [10$^6$ m/s]')
ax[0].set_title('electrons at t = 0: uniform in x, Maxwellian at %.0f K (%.2f eV)' % (Te0, te_ev), fontsize=10)
for k in (-1, 1):
    ax[0].axhline(k * vte * 1e-6, color='k', ls='--', lw=0.8)
ax[0].text(0.2, vte * 1e-6 * 1.1, r'$\pm v_{th,e}=\sqrt{kT_e/m_e}$ = %.2e m/s' % vte, fontsize=9)
sub = np.random.default_rng(1).choice(len(io), size=min(20000, len(io)), replace=False)
ax[1].plot(io[sub, 1] * 1e2, io[sub, 2] * 1e-3, ',', color='C3', alpha=0.6)
ax[1].set_xlabel('x [cm]'); ax[1].set_ylabel('v$_x$ [10$^3$ m/s]')
ax[1].set_title('He$^+$ at t = 0: Maxwellian at %.0f K (%.3f eV)' % (Ti0, KB * Ti0 / QE), fontsize=10)
for k in (-1, 1):
    ax[1].axhline(k * vti * 1e-3, color='k', ls='--', lw=0.8)
ax[1].text(0.2, vti * 1e-3 * 1.1, r'$\pm v_{th,i}$ = %.0f m/s' % vti, fontsize=9)
for a in ax:
    a.axvline(0, color='gray'); a.axvline(L * 1e2, color='gray')
fig.suptitle('%d macro-particles per species, weight w = %.3e m$^{-2}$ each, %d per cell' % (len(el), w, int(info['nppc'])), fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(out, 'init_phase_space.png'), dpi=130); plt.close(fig)

# ---- 2. velocity and energy distributions ----
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
for a, v, vth, lab, col in [(ax[0], el[:, 2], vte, 'electrons', 'C0'), (ax[1], io[:, 2], vti, 'He$^+$', 'C3')]:
    h, edges = np.histogram(v, bins=80, density=True)
    cen = 0.5 * (edges[1:] + edges[:-1])
    a.semilogy(cen / vth, h * vth, '.', color=col, label='loaded particles')
    vv = np.linspace(-5, 5, 400)
    a.semilogy(vv, np.exp(-vv**2 / 2) / np.sqrt(2 * np.pi), 'k-', lw=1, label='Maxwellian')
    a.set_ylim(1e-5, 1); a.set_xlabel(r'$v_x / v_{th}$'); a.set_ylabel(r'$f(v_x)\,v_{th}$')
    a.set_title(lab + ': v$_x$ distribution'); a.legend()
e_ev = 0.5 * ME * (el[:, 2]**2 + el[:, 3]**2 + el[:, 4]**2) / QE
h, edges = np.histogram(e_ev, bins=80, range=(0, 25), density=True)
cen = 0.5 * (edges[1:] + edges[:-1])
ax[2].semilogy(cen, h, '.', color='C0', label='loaded electrons')
ee = np.linspace(0.01, 25, 400)
ax[2].semilogy(ee, 2 * np.sqrt(ee / np.pi) * te_ev**-1.5 * np.exp(-ee / te_ev), 'k-', lw=1,
               label=r'$2\sqrt{\epsilon/\pi}\,T_e^{-3/2}e^{-\epsilon/T_e}$')
for eth, name in [(19.82, 'exc. 19.82'), (20.61, 'exc. 20.61'), (24.59, 'ioniz. 24.59')]:
    ax[2].axvline(eth, color='gray', ls=':', lw=0.8)
    ax[2].text(eth, 3e-1, name, rotation=90, fontsize=7, va='top', ha='right')
ax[2].set_ylim(1e-7, 1); ax[2].set_xlabel('electron energy [eV]'); ax[2].set_ylabel('f($\\epsilon$) [eV$^{-1}$]')
ax[2].set_title('electron energy distribution, <$\\epsilon$> = %.2f eV = 3/2 T$_e$' % e_ev.mean()); ax[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(out, 'init_velocities.png'), dpi=130); plt.close(fig)

# ---- 3. grid quantities ----
fig, ax = plt.subplots(2, 2, figsize=(12, 7.5))
a = ax[0, 0]
a.plot(xg * 1e2, grid[:, 4], label='n$_e$ (deposited)')
a.plot(xg * 1e2, grid[:, 5], label='n$_i$ (deposited)')
a.axhline(n0, color='k', ls='--', lw=0.8, label='n$_0$')
sig = n0 / np.sqrt(info['nppc'])
a.fill_between(xg * 1e2, n0 - sig, n0 + sig, color='gray', alpha=0.2, label=r'$\pm n_0/\sqrt{n_{ppc}}$ shot noise')
a.set_xlabel('x [cm]'); a.set_ylabel('m$^{-3}$'); a.legend(fontsize=8); a.set_title('densities at t = 0 (CIC deposit of the random loading)')
a = ax[0, 1]
a.plot(xg * 1e2, grid[:, 1] / QE, 'k')
a.set_xlabel('x [cm]'); a.set_ylabel(r'$\rho/e$ [m$^{-3}$]'); a.set_title(r'charge density: $\rho = e(n_i - n_e)$, pure noise, rms = %.2e m$^{-3}$' % (grid[1:-1, 1].std() / QE), fontsize=10)
a = ax[1, 0]
a.plot(xg * 1e2, grid[:, 2], 'C2')
a.set_xlabel('x [cm]'); a.set_ylabel(r'$\phi$ [V]'); a.set_title(r'potential at t = 0: $\phi(0)=V(0)=0$, $\phi(L)=0$, noise rms = %.3f V' % grid[:, 2].std(), fontsize=10)
a = ax[1, 1]
a.plot(xg * 1e2, grid[:, 3], 'C1')
a.set_xlabel('x [cm]'); a.set_ylabel('E [V/m]'); a.set_title('field at t = 0: rms = %.1f V/m  (sheath fields later: ~10$^4$-10$^5$ V/m)' % grid[:, 3].std(), fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(out, 'init_grid.png'), dpi=130); plt.close(fig)

# ---- 4. drive and scales ----
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
T = 1 / freq
t = np.linspace(0, T, 1000)
ax[0].plot(t / T, V0 * np.sin(2 * np.pi * freq * t), 'k')
spp = int(info['steps_per_period'])
tk = np.arange(0, spp + 1) * dt
ax[0].plot(tk / T, V0 * np.sin(2 * np.pi * freq * tk), '.', ms=2, color='C3', label='%d time steps per period, dt = %.3e s' % (spp, dt))
nph = int(info['n_phase'])
for k in range(0, nph + 1, nph // 10):
    ax[0].axvline(k / nph, color='gray', lw=0.5, ls=':')
ax[0].set_xlabel('t / T   (T = 1/f = %.2f ns)' % (T * 1e9)); ax[0].set_ylabel(r'$\phi(0,t)$ [V]')
ax[0].set_title(r'drive: $\phi(0,t) = V_0 \sin(2\pi f t)$, V$_0$ = %.0f V, f = %.2f MHz; grounded at x = L' % (V0, freq * 1e-6), fontsize=10)
ax[0].legend(fontsize=8)
a = ax[1]
lengths = [('L (gap)', L), ('dx', dx), (r'$\lambda_D$ (t=0)', lam_d), (r'$\lambda_D$ (steady, n=1.4e14, T=2.5 eV)',
            np.sqrt(EPS0 * 2.5 * QE / (1.4e14 * QE**2)))]
times = [('T = 1/f', T), ('dt', dt), (r'$1/\omega_{pe}$ (t=0)', 1 / wpe), (r'ion transit L/2c$_s$', 0.5 * L / np.sqrt(KB * Te0 / MION)),
         (r'$\omega_{pe}^{-1}$ steady', 1 / np.sqrt(1.4e14 * QE**2 / (EPS0 * ME)))]
if ngas > 0:
    sig_tot = 5e-20
    lengths.append((r'e mean free path $1/(n_g\sigma)$, $\sigma$=5e-20', 1 / (ngas * sig_tot)))
    times.append((r'1/$\nu_{max,e}$ (8.65e7 s$^{-1}$)', 1 / 8.65e7))
names = [n for n, _ in lengths] + [n for n, _ in times]
vals = [v for _, v in lengths] + [v for _, v in times]
cols = ['C0'] * len(lengths) + ['C3'] * len(times)
y = np.arange(len(names))
a.barh(y, vals, color=cols)
a.set_xscale('log'); a.set_yticks(y); a.set_yticklabels(names, fontsize=8); a.invert_yaxis()
for yi, v in zip(y, vals):
    a.text(v, yi, ' %.2e' % v, va='center', fontsize=8)
a.set_title(r'scales: lengths [m] (blue), times [s] (red); $\lambda_D$/dx = %.2f, $\omega_{pe}$dt = %.3f at t = 0' % (lam_d / dx, wpe * dt), fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(out, 'init_drive.png'), dpi=130); plt.close(fig)

# ---- 5. cross sections and collision frequencies ----
here = os.path.dirname(os.path.abspath(__file__))
dhe = os.path.join(here, '..', 'data', 'he')
if os.path.isdir(dhe):
    ng = ngas if ngas > 0 else 9.64e20
    tabs = {}
    for name in ['electron_scattering', 'excitation_1', 'excitation_2', 'ionization', 'ion_scattering', 'ion_back_scatter']:
        tabs[name] = np.loadtxt(os.path.join(dhe, name + '.dat'))
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    a = ax[0]
    for name, lab in [('electron_scattering', 'e elastic'), ('excitation_1', 'e exc. 19.82 eV'), ('excitation_2', 'e exc. 20.61 eV'), ('ionization', 'e ionization 24.59 eV')]:
        t_ = tabs[name]; a.loglog(t_[:, 0], t_[:, 1], label=lab)
    a.set_xlim(1e-2, 1e3); a.set_ylim(1e-24, 1e-18); a.set_xlabel('electron energy [eV]'); a.set_ylabel(r'$\sigma$ [m$^2$]'); a.legend(fontsize=8); a.set_title('e + He cross sections (Biagi 7.1, via WarpX data)')
    a = ax[1]
    for name, lab in [('ion_scattering', 'He$^+$ isotropic'), ('ion_back_scatter', 'He$^+$ backward')]:
        t_ = tabs[name]; a.loglog(t_[:, 0], t_[:, 1], label=lab)
    a.set_xlim(1e-3, 1e3); a.set_xlabel(r'$\epsilon_{rel} = \frac{1}{2}\mu g^2$ [eV], $\mu$ = M/2'); a.set_ylabel(r'$\sigma$ [m$^2$]'); a.legend(fontsize=8); a.set_title('He$^+$ + He cross sections (Phelps)')
    a = ax[2]
    E = tabs['electron_scattering'][:, 0]
    E[0] = 1e-3
    v = np.sqrt(2 * E * QE / ME)
    tot = np.zeros_like(E)
    for name in ['electron_scattering', 'excitation_1', 'excitation_2', 'ionization']:
        t_ = tabs[name]; s_ = np.interp(E, t_[:, 0], t_[:, 1], left=0, right=t_[-1, 1]); tot += s_
        a.loglog(E, ng * s_ * v, lw=0.8, label=name.replace('_', ' '))
    nu = ng * tot * v
    i = np.argmax(nu)
    a.loglog(E, nu, 'k', lw=2, label=r'total $\nu = n_g \sigma_{tot} v$')
    a.axhline(nu[i], color='r', ls='--', label=r'$\nu_{max}$ = %.3e s$^{-1}$ at %.1f eV' % (nu[i], E[i]))
    a.set_xlim(1e-2, 1e3); a.set_ylim(1e5, 2e8); a.set_xlabel('electron energy [eV]'); a.set_ylabel(r'$\nu$ [s$^{-1}$]')
    a.set_title(r'electron collision frequencies at n$_g$ = %.2e m$^{-3}$; P = 1 - e$^{-\nu_{max} dt}$ = %.2f %%' % (ng, 100 * (1 - np.exp(-nu[i] * dt))), fontsize=10)
    a.legend(fontsize=7, loc='lower right')
    fig.tight_layout(); fig.savefig(os.path.join(out, 'init_xsec.png'), dpi=130); plt.close(fig)

print('figures written to', out)
print('summary: n0 = %.3e m^-3, Te0 = %.3f eV, lambda_D = %.3e m (lambda_D/dx = %.2f), wpe = %.3e rad/s (wpe dt = %.3f), v_th,e dt/dx = %.3f'
      % (n0, te_ev, lam_d, lam_d / dx, wpe, wpe * dt, vte * dt / dx))
print('loaded: <vx^2>^(1/2)/v_th,e = %.4f (electrons), %.4f (ions); mean electron energy %.4f eV (3/2 Te = %.4f eV)'
      % (np.sqrt((el[:, 2]**2).mean()) / vte, np.sqrt((io[:, 2]**2).mean()) / vti, e_ev.mean(), 1.5 * te_ev))
print('grid noise at t = 0: n_e rms/n0 = %.4f (1/sqrt(nppc) = %.4f), phi rms = %.4f V, E rms = %.2f V/m'
      % (grid[1:-1, 4].std() / n0, 1 / np.sqrt(info['nppc']), grid[:, 2].std(), grid[:, 3].std()))
