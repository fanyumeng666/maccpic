# Data files for ccp1d

## he/ — electron–He and He⁺–He cross sections

Copied unchanged on 2026-11-01 from the WarpX data repository
`https://github.com/BLAST-WarpX/warpx-data`, directory `MCC_cross_sections/He/`
(licence: Community Data License Agreement – Permissive 2.0, see
`he/LICENSE_CDLA-Permissive-2.0.txt`). These are the tables used by the WarpX
capacitive-discharge example that reproduces the Turner et al. (2013) benchmark.

Original sources, as stated in the warpx-data README:

| file | process | source |
| --- | --- | --- |
| electron_scattering.dat | e + He elastic (isotropic) | S. F. Biagi, MAGBOLTZ 7.1 (June 2004), via LXCat |
| excitation_1.dat | e + He excitation, threshold 19.82 eV | Biagi 7.1 |
| excitation_2.dat | e + He excitation, threshold 20.61 eV | Biagi 7.1 |
| ionization.dat | e + He → 2e + He⁺, threshold 24.59 eV | Biagi 7.1 |
| ion_scattering.dat | He⁺ + He isotropic elastic | A. V. Phelps database (LXCat) |
| ion_back_scatter.dat | He⁺ + He backward scattering | A. V. Phelps database (LXCat) |
| charge_exchange.dat | He⁺ + He charge exchange, σ_ct = σ_bs + 0.5 σ_iso | derived; **not used** in the Turner benchmark |

Format: two columns, `energy [eV]   cross section [m^2]`, 0.01 eV spacing, up to
749.99 eV. Beyond the last tabulated energy use the last value (Turner's
convention). Energy variable: electron kinetic energy for the electron files;
**centre-of-mass (relative) energy ε_rel = ½ μ g², μ = M_He/2** for the ion
files (WarpX looks all MCC cross sections up at the reduced-mass collision
energy; Turner 2013 defines the ion processes in the centre-of-mass frame).

Thresholds in the files: excitation 19.82 and 20.61 eV, ionization 24.59 eV.

## turner2013_case1_ni.dat

Time-averaged ion density of Turner et al., Phys. Plasmas 20, 013507 (2013),
benchmark case 1 (He, 450 V, 13.56 MHz, gap 6.7 cm, n_gas 9.64e20 m⁻³), on the
129 grid nodes, as tabulated in the WarpX example script
`Examples/Physics_applications/capacitive_discharge/analysis_1d.py` (taken from
the paper's supplementary data). Columns `x [m]   n_i [m^-3]`. Centre value
1.4046e14 m⁻³, edge value 3.6e13 m⁻³. The five codes in the paper agree with
each other to better than 1 %.
