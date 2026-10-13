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
