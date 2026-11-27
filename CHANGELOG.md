# Changelog

All notable changes to libnest are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/) with the 0.x caveat that the public API
is still settling.

Releases 0.1.0 and 0.1.1 predate this file; see the git history for their contents.

## [Unreleased]

Work in progress on the `edf-refactor` branch (see `TODO_edf.md`): a single
implementation of the Skyrme formulas behind `libnest.edf`, selectable
parametrizations, and fixes for the unit errors listed under "Known issues" below.

## [0.2.0] — 2026-11-27

The mean field `U_q` now follows W-BSK: it depends on tau and the currents, and
optionally on the density gradients. Its values change (symmetric matter at
rho_n = rho_p = 0.08 fm^-3: -59.8 MeV; 0.1.2 gave -18.4 MeV). Existing calls
`U_q(rho_n, rho_p, q)` still work; the new arguments are optional.

### Fixed — results differ from 0.1.2

- **`bsk.U_q` and `bskg.U_q`: wrong t3 term.** The mean field is the derivative of
  the t3 energy t3/12 rho^alpha [(1 + x3/2) rho^2 - (x3 + 1/2)(rho_n^2 + rho_p^2)],
  i.e. t3/12 rho^(alpha-1) [(1 + x3/2)(alpha+2) rho^2 - (x3 + 1/2)(2 rho rho_q
  + alpha (rho_q^2 + rho_q'^2))]. The code had `(0.5+0.5*X3)` instead of
  `(1+0.5*X3)` and a product `2*rho*rho_q*ALPHA*rho_q'**2` instead of the sum, so the
  t3 repulsion was too strong: at rho_n = rho_p = 0.08 fm^-3 U_n was -18.4 instead of
  -65.8 MeV, and in pure neutron matter at 0.16 fm^-3 +198.5 instead of -37.1 MeV.
  `U_q` now equals d(`epsilon_rho_np`)/d(rho_q) (new test `TestMeanField`). It still
  contains only the t0, t3 terms: no d(epsilon_tau)/d(rho_q), gradient or Coulomb
  terms. Affects `plots.plot_U_q` and `real_data_plots`.
- **`bsk.U_q` and `bskg.U_q`: no tau terms.** `U_q` is now `U_rho_q + U_tau_q`, as
  `U_rho_q + U_tau_q` of W-BSK (`hpc-engine/bsk/bsk_functional.c`, nuclear matter):
  d(`epsilon_rho_np` + `epsilon_tau_np`)/d(rho_q) at fixed tau, j. The tau part comes from
  the density dependence of C^tau_0, C^tau_1 (t4, t5) and from C^tau_1 (rho_n - rho_p)(tau_n -
  tau_p): U_tau_q = C0_tau (tau_n + tau_p) + C1_tau (tau_q - tau_q') + dC0_tau [rho (tau_n +
  tau_p) - (j_n + j_p)^2] + dC1_tau [(rho_n - rho_p)(tau_n - tau_p) - (j_n - j_p)^2].
  New signature `U_q(rho_n, rho_p, q, tau_n=None, tau_p=None, jsum2=0, jdiff2=0)`; without
  tau the uniform Fermi gas `definitions.rho2tau(rho_q)` is used. Symmetric matter at
  rho_n = rho_p = 0.08 fm^-3: U_q = -59.8 MeV (was -65.8), pure neutron matter at 0.16
  fm^-3: U_n = -24.5 MeV (was -37.1). Checked against the C code of W-BSK: equal to 1e-15
  with DENSEPSILON = 0 on both sides (W-BSK has 1e-8 inside the powers, libnest 1e-12:
  ~1e-5 MeV apart at these densities). Still not included (in `V_q` of W-BSK): pairing
  rearrangement, Coulomb. The gradient terms are optional, see below.
- **`bsk.epsilon_delta_rho_np`: wrong cross term.** grad rho_n . grad rho_p was computed as
  `0.5*(p² - n² - p²)` = -|grad rho_n|²/2 instead of `0.5*(|grad rho|² - n² - p²)`; this
  changes the t4 (beta) part. Now equal to `e_LaplaceRho_np` of W-BSK (1e-16 with
  DENSEPSILON = 0). `epsilon_np` passed |grad rho_n|² + |grad rho_p|² as |grad rho|²; it now
  uses (grad rho_n + grad rho_p)². The same in `bskg`.

### Added

- `bsk.U_rho_q`, `bsk.U_tau_q` (the two parts of `U_q`) and `bsk.dC0_rho`, `dC1_rho`,
  `dC0_tau`, `dC1_tau` (d/d rho of the coupling constants); the same in `bskg`.
- `bsk.U_LaplaceRho_q(rho_n, rho_p, q, rho_grad_n_square, rho_grad_p_square,
  rho_grad_square, rho_lap_n, rho_lap_p)`: the gradient part of the mean field, the functional
  derivative of int `epsilon_delta_rho_np` over rho_q, ported from `U_LaplaceRho_n/p` of W-BSK
  (equal to machine precision with DENSEPSILON = 0; checked against a numerical functional
  derivative on a grid). `U_q` takes the same five gradient arguments as optional keywords and
  then adds this term (all five or none). In static Pb208 / Ca40 BSk31 solutions from MOCCa
  the weighted rms of V_q - U_q drops from 3.7 / 6.6 MeV (rho, tau only) to 1.1 / 0.27 MeV
  with the gradient term; the rest is mostly the spin-orbit term, which libnest does not
  have. The same in `bskg`. libnest does not compute derivatives: the gradients and
  Laplacians are inputs, computed on the user's grid (e.g. by FFT).

### Known issues

- `U_q` does not contain the pairing rearrangement term, the spin-orbit term, or the
  Coulomb potential for protons.
- `definitions.mu_q` and `bsk.v_pi` still have the unit error listed under 0.1.2, and
  so do `I` and `epsilon_pi_np`, which are built on them.

## [0.1.2] — 2026-10-13

### Changed — results differ from 0.1.1

- **`definitions.rho2tau` was low by a factor of pi^(2/3) ~ 2.15.** The kinetic
  density of a uniform Fermi gas is tau = (3/5) k_F^2 rho = (3/5)(3 pi^2)^(2/3)
  rho^(5/3); the code used `(3*pi)^(2/3)`. At rho = 0.08 fm^-3 it now returns
  0.0853 instead of 0.0398 fm^-5. Anything built on `rho2tau` changes accordingly,
  including the `tau_vs_kf` / `tau_vs_rho` figures.
- **Reference pairing fields** now evaluate the symmetric-matter gap at half the
  total density where the interpolation requires it, and are clamped to
  non-negative values, so an asymmetric-matter interpolation can no longer return a
  negative gap.

### Added

- Asymmetric-matter pairing interpolation schemes, with unit tests:
  `neutron_ref_pairing_field_eq2` / `_eq3` / `_eq6` and the matching
  `proton_ref_pairing_field_*` (Chamel 2009; BSkG3 Eq. 3; BSkG4 Eq. 6).
  `libnest.bsk` still defaults to the Chamel 2009 scheme.
- `libnest/bskg.py` — a placeholder for the BSkG family. It is a copy of `bsk.py`
  kept for in-progress work and is **not** a separate functional yet; do not depend
  on it.
- A strict documentation build in CI (`-W --keep-going -n`), so any Sphinx warning
  or unresolved cross-reference fails the build.

### Fixed — documentation and conventions

- `rho2kf` / `kf2rho`: the stated convention was wrong. rho is the density of one
  isospin with **both** spin components, which is how the code is used throughout.
  `kf2rho` also had Args and Returns swapped.
- Effective masses (`effMn`, `effMp`, `isoscalarM`, `isovectorM`, `B_q`) are
  documented as the dimensionless ratio M*/M; the docstrings claimed MeV.
- `units`: nucleon-mass table is 10^-27 kg (was 10^27); `KtoMev` had the conversion
  inverted; `gcm3tofm3` described it backwards; `RHOND` is 0.00024 fm^-3 (was
  0.00042); `RHOSAT` is stated as the crust-core density 0.179 fm^-3, not the
  nuclear saturation density n_0; `hbar22M0` documented.
- The "Bulk Neutron Properties" table on the Inner Crust page renders again (a CSV
  quoting error had hidden it); the Landau-velocity figures are no longer blank;
  five bibliography entries cited in `physics.rst` but missing from the `.bib` were
  added; duplicate-citation warnings removed in favour of one global bibliography.
- Development TODO boxes and the stray `topsecret` page are no longer published.
- The installation page describes `pip install libnest` instead of a `git clone`
  with a `wdata` requirement; `instalation.rst` is now `installation.rst`.
- Hard-coded absolute paths removed from the examples and the figure scripts; the
  figure pipeline rebuilds when the library changes, not only when a script does.
- About 50 broken `:func:` cross-references and wrong type names in docstrings.

### Packaging

- `project.license` is an SPDX string and the deprecated license classifier is
  gone (setuptools >= 77); this silences a deprecation that becomes an error in
  February 2027.

### Known issues

Unit errors in the chemical potential and the pairing strength, present in this
release and fixed on `edf-refactor` (see `TODO_edf.md` section 5):

- `definitions.mu_q` treats the dimensionless ratio M*/M as a mass in MeV, so it
  returns about 31000 MeV where a few tens of MeV are expected.
- `bsk.v_pi` has the same error and should use B_q^(3/2).
- The t3 term of `bsk.U_q` is dimensionally inconsistent, and the gradient terms of
  `epsilon_np` / `epsilon_delta_rho_np` drop a cross term.

`I`, `epsilon_pi_np` and `epsilon_np` are built on these and are affected too.
