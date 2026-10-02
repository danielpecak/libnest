# TODO_doc — Documentation & website roadmap

Goal: the published site (`danielpecak.github.io/libnest`, built by Sphinx from `docs/`)
should be **correct, complete where it is visible, and impossible to silently rot**.
A student should land on it, install with `pip install libnest`, run a 5-line example,
and find every public function documented with units and conventions.

Companion plan: [`TODO_edf.md`](TODO_edf.md) (multi-functional refactor). The API-reference
work below (D3) is scheduled **after** EDF phase E2, because the BSk docstrings move into a
new class there. Everything else here can start now.

Priorities: **D0** = wrong/broken on the live site · **D1** = build & deploy infrastructure
(make breakage visible) · **D2** = information architecture · **D3** = API reference quality ·
**D4** = new content · **D5** = polish.

---

## How this analysis was done

- Built the docs locally (Sphinx 9.1, scratch venv, figures regenerated from the scripts):
  **38 lines of warnings/errors** in a normal build, **~120** with nitpicky mode (`-n`),
  which catches broken `:func:` cross-references.
- Viewed every `.rst` page, the module docstrings rendered by `autodoc`, every generated figure,
  `conf.py`, both Makefiles and `.github/workflows/documentation.yml`.
- Spot-checked the numbers behind figures and formulas against the code.

---

## D0 — Wrong or broken on the live site (fix first; small diffs)

_Status (2026-10-02): done except "Collisions Initial State", which needs a corrected
description from you. Result: the normal Sphinx build has **0 warnings** (was 38 lines);
nitpicky mode still reports ~85, all docstring cross-reference issues scheduled for D3._

Done in this pass, beyond the items' original wording:
- `rho2tau` fixed to (3π²)^(2/3); `test_rho2tau_basic` now checks τ = 3/5 k_F² ρ (+ array test).
- `units`: masses 10⁻²⁷ kg; `KtoMev` formula; `gcm3tofm3` described the conversion backwards;
  `RHOND` table said 0.00042 fm⁻³ (actual 0.00024); `RHOSAT` 0.179 fm⁻³ also in README;
  `hbar22M0` documented. `CLAUDE.md` convention for `rho2kf` corrected.
- Bibliography: the 5 missing entries added (metadata checked against Crossref); per-page
  `.. bibliography::` removed from `bsk.rst`/`physics.rst` → one global bibliography page.
- Pages removed: `topsecret.rst`, `plotting.rst` (all sections empty); `instalation.rst` →
  `installation.rst`; `help.rst` → "Citing and Contributing"; `pasta.rst` marked as placeholder.
- Hard-coded paths removed: `docs/_static/slices.py` → `examples/wdata_slices.py` (files from
  the command line); `examples/tools_example.py` takes the file as an argument. No absolute
  paths remain in the repo (`CLAUDE.md` updated accordingly).

### Physics statements that are wrong (need your sign-off; touch code + tests)
- [x] **`rho2tau` is off by a factor π^(2/3) ≈ 2.15.** Code and docstring use
      `0.6*(3π)^(2/3) ρ^(5/3)` (`libnest/definitions.py`, `rho2tau`); the uniform Fermi-gas
      result is τ = (3/5) k_F² ρ = (3/5)(3π²)^(2/3) ρ^(5/3). At ρ = 0.08 fm⁻³:
      libnest 0.0398, textbook 0.0853 fm⁻⁵. The figures `tau_vs_kf.png` / `tau_vs_rho.png`
      show the wrong curve. `tests/test_definitions.py::test_rho2tau_basic` re-implements the
      same formula, so it locks the bug in. Replace that test with the invariant τ = 3/5 k_F² ρ.
- [x] **`rho2kf` / `kf2rho` docstrings state the wrong convention.** They say ρ is "the
      density of one isospin and one spin component only". The formula k_F = (3π²ρ)^(1/3)
      is for **one isospin, both spin components** — which is also how the code is used
      everywhere (benchmark: `rho2kf(0.08) = 1.33` in SNM; the C code does the same). The same
      wrong sentence is in `CLAUDE.md` ("Physics conventions"). `kf2rho` also has Args/Returns
      swapped (says the input is a density and the output a wavevector).
- [x] **`units` module docstring:** nucleon-mass table says `·10^27 kg` (should be `10^-27`);
      `KtoMev` shows E ≈ 1.16·10¹⁰ · T (it is T / 1.16·10¹⁰); the neutron-star table gives
      RHOSAT = 0.18 fm⁻³ while README lists 0.16 — state clearly that `RHOSAT` is
      3·10¹⁴ g/cm³ ≈ 0.179 fm⁻³ (crust–core), not the nuclear saturation density n₀ ≈ 0.16.
      `hbar22M0 = 20.72` is undocumented and duplicates `HBAR2M_n`.
- [x] **Effective-mass docstrings claim `[MeV]`** but `effMn`, `effMp`, `isoscalarM`,
      `isovectorM` return the dimensionless ratio M*/M. `B_q` says "Returns: effective mass of a
      proton". (`mu_q`/`v_pi` treat that ratio as a mass in MeV — real bug, tracked in
      `TODO_edf.md`.)

### Visibly broken pages / figures
- [x] **First table on "Inner Crust" is not rendered at all** (docutils CSV error,
      `docs/inner-crust.rst:17`): the header has `":math:`\\epsilon_{\\mathrm{F}}`" [MeV]`
      (quote closed too early) and `[\MeV]`. Fix the quoting; the whole "Bulk Neutron
      Properties" table then reappears.
- [x] **Landau-velocity figures are empty.** `vLandau` returns a fraction of c (max ≈ 0.013);
      `plot_vlandau_vs_kf.py` / `_vs_rho.py` set `ylim([0, 2.7])` and label the axis "% c",
      so the curve is a flat line at the bottom. Multiply by 100 (and keep "% c") or fix the
      limits. Also start `kf` at > 0 (divide-by-zero warnings during the build).
- [x] **TODO boxes are published.** `conf.py` has `todo_include_todos = True`, so every module
      page shows "Describe properties / Give references" boxes, and the "Other" page lists 15
      of them. Set it to `False` for the public build (optionally enable via an env var for
      local builds) and move the real TODOs into this file.
- [x] **`topsecret.rst` is built and published** (`topsecret.html` is reachable) although it
      is not in any toctree. Delete it, or add it to `exclude_patterns`.
- [x] **Bibliography:** 5 keys cited in `physics.rst` are missing from `bibtexNS.bib`
      (`chamel2008further`, `chamel2010spin`, `goriely2013further`, `goriely2013hartree`,
      `goriely2015further`), so the "History of BSk" list renders unresolved citations.
      13 "duplicate citation" warnings come from having a `.. bibliography::` in `bsk.rst`,
      `physics.rst` *and* `bibliography.rst`. Choose one model: a global bibliography page
      + `:cite:` everywhere, **or** per-page `:footcite:` + `.. footbibliography::`.
- [ ] **"Collisions Initial State"** (`inner-crust.rst`) shows a red
      `.. error:: The system description is not accurate.` box on the public page. Fix the
      description or remove the section until it is correct.

### Outdated or placeholder content
- [x] **Installation page is wrong.** It says the library *requires* `wdata` and must be
      installed with `git clone`. Reality: `pip install libnest`; runtime deps
      NumPy/SciPy/Matplotlib/pandas; `wdata` only for the WDATA tools. Rename
      `instalation.rst` → `installation.rst` (typo in the URL as well).
- [x] **Tutorial is a placeholder** ("1 2 3") whose only content is a `literalinclude` of
      `_static/slices.py` — a script with a hard-coded `/home/pecak/sshfs/...` path that needs
      `wdata` and private data. Replace (see D2).
- [x] **Empty headings on public pages:** `help.rst` (titled "Other": empty Citations,
      Development, Contributing), `plotting.rst` (every section empty), `physics.rst`
      ("Pairing P: Give some references.", empty "Pasta phase"), `pasta.rst` /
      `io.rst` / `tools.rst` (empty "Minkowski Functionals", "Regular TXT files", "WDATA").
      Rule: no empty heading ships — fill it or remove it.
- [x] **Stale pointers outside `docs/`:** `main.py` tells users to see
      `examples/legacy_tests.py` (deleted in P3); `examples/tools_example.py` has a hard-coded
      `/media/data/...` path and needs `wdata` — label it as such or rework it.
- [x] **Typos visible on the site:** "genealized", "within fully dynamical approach",
      "desnity" (×2, `units`), "neuton", "pairng", "voume", "excesive", "occure",
      "superlufid", "nergy", "BskG4"→"BSkG4".

---

## D1 — Build & deploy infrastructure (make breakage visible)

- [ ] **CI installs the wrong things.** In `documentation.yml`,
      `pip install sphinx>=7.0 sphinx_rtd_theme ...` is unquoted: the shell treats `>=7.0` as
      an output redirect (creates a file named `=7.0`) and the version floor is ignored.
      Replace both `pip install` lines with `pip install -e ".[docs]"` — this also installs
      libnest itself, so autodoc and the figure scripts import the installed package.
- [ ] **Fail the build on warnings.** After D0, switch to
      `sphinx-build -W --keep-going -n` (warnings are errors, nitpicky cross-refs). This is
      the single change that stops the site from rotting again.
- [ ] **`conf.py` cleanup:**
      - add `intersphinx_mapping` for python/numpy/scipy/matplotlib/pandas (the extension is
        loaded but has no mapping — `np.gradient` refs are unresolved);
      - add `sphinx.ext.viewcode` (source links) and `sphinx.ext.doctest` (see D3);
      - `autodoc_member_order = "bysource"` (functions currently appear alphabetically, which
        scatters related functions);
      - remove `templates_path = ['_templates']` (directory doesn't exist) or create it;
      - replace the Texinfo placeholder `'One line description of project.'`;
      - align `needs_sphinx` with the `>=7.0` floor.
- [ ] **Figure pipeline (`docs/source/Makefile`):**
      - figures depend only on their `plot_*.py`, not on `libnest/*.py`, so after a physics fix
        `make` keeps the stale PNGs → add the library sources as prerequisites;
      - `clean` uses `rm` without `-f` (fails if a figure is missing); `mkdir -p` runs *after*
        the figures are written;
      - drop `sys.path.insert(0, '../../')` and `from libnest.bsk import *` from the scripts
        once CI installs the package;
      - add a tiny shared style helper (`plots/_style.py`: figure size, dpi, axis labels with
        units, `plt.close()`), so all figures look like one set.
      - _Option to evaluate:_ matplotlib's built-in `.. plot::` directive
        (`matplotlib.sphinxext.plot_directive`) keeps the "figures come from code" rule, puts
        the code next to the figure in the `.rst`, and removes the Makefile layer entirely.
- [ ] **Repo hygiene:** delete the stray `_sources/physics.rst.txt` at the repo root (left
      over from a 2024 `sphinx`↔`main` merge), the empty `docs/_static/workaround.txt` and
      `docs/source/plots/TODO.txt` (`slices.py` already moved to `examples/` in D0);
      move `docs/DOCUMENTATION_REVIEW_2025-11-19.md` to `docs/history/` (its open items
      are absorbed here); update `docs/README_DOCS.md` (install via `.[docs]`, figure pipeline,
      `-W`).
- [ ] **One docs URL everywhere.** `pyproject.toml` points `Documentation` to
      `libnest.readthedocs.io` (dead); README uses GitHub Pages. Pick GitHub Pages and confirm
      Settings → Pages serves the `sphinx` branch.
- [ ] **Link check:** add `make linkcheck` as a non-blocking CI step (or scheduled). README's
      "Related projects → SkyNET" link looks unrelated/wrong.

---

## D2 — Information architecture (what the reader sees)

Today the left sidebar is "Documentation / Modules / Examples", with module pages that are
mostly raw `automodule` dumps. Proposed structure:

```
Home (index)        short pitch · pip install · 5-line example · one figure · links
Getting started
  Installation
  Quickstart          README "Quick Start", executed as a doctest
  Units & conventions fm⁻³ / MeV / fm; per-species density (both spins); DENSEPSILON;
                      scalar + array inputs; ħ²/2m per species  ← most useful single page
User guide
  Uniform matter: EOS, pressure, speed of sound
  Pairing in asymmetric matter (interpolation schemes, comparison figure)
  Effective masses and mean fields
  Analyze your own data   (examples/analyze_data.py + its generated figure)
  Working with WDATA      (optional dependency)
Physics background
  Neutron-star structure · Skyrme & BSk functionals · Pairing models
  Inner-crust data tables · Bibliography
API reference
  one page per public module, autosummary tables grouped by topic
  "Research scripts (not portable)": real_data_plots, delta_and_temperature, mass_center, nucleus
Development
  Contributing · Testing philosophy · Building docs & figures · Release process
  Changelog · How to cite (CITATION.cff) · Contributors & funding
```

- [ ] Rewrite `index.rst` as a landing page; move Contributors/Funding/LUMI to an "About" page.
- [ ] Create **Units & conventions** (pull from `units` docstring, `CLAUDE.md` "Physics
      conventions", and the corrected `rho2kf` convention).
- [ ] Replace `tutorial.rst` / `plotting.rst` with the User-guide pages above; every page runs
      only in-repo code and produces in-repo figures.
- [ ] Split `plots.rst`: `libnest.plots` stays in the API reference; `libnest.real_data_plots`
      (hard-coded data paths, needs private files) moves under "Research scripts" with a
      warning box, or is dropped from the public docs.
- [ ] Replace `help.rst` with Development pages (content exists in README "Contributing" and
      in `CLAUDE.md` "Testing philosophy").

---

## D3 — API reference quality (after EDF phase E2)

Docstrings are the bulk of the site. Do this pass **once**, on the new `SkyrmeFunctional`
class from `TODO_edf.md`, not on `bsk.py` functions that are about to become thin wrappers.

- [ ] **Topic-grouped API pages generated with `autosummary`** instead of the hand-maintained
      lists in `bsk.rst`, which already went stale: they reference non-existent
      `E_minigap_n`, `derivative_pressure_rho_n`, `derivative_epsilon_rho_n`, list
      `definitions` functions (`eF_n`, `mu_q`) under a heading called "Move", and omit the
      `_eq2/_eq3/_eq6` pairing schemes.
- [ ] **Fix ~20 broken cross-references in `plots.py`** (`:func:\`energy_per_nucleon\`` needs a
      leading dot or the full path) and in `real_data_plots.py`; replace type names `string`
      → `str`, `numpy` → `numpy.ndarray`.
- [ ] **Remove the duplicate page titles:** module docstrings start with
      `Module: X\n=====` and the `.rst` files repeat the same title, so it appears twice in the
      sidebar/TOC. Keep the title in the `.rst` only.
- [ ] **Don't document what doesn't work:** `bsk.testMe` (debug helper), and in `plots.py`
      `plot_epsilon`, `plot_epsilon_tau`, `plot_epsilon_delta`, `plot_epsilon_rho_np`,
      `plot_epsilon_tau_np`, `plot_epsilon_delta_rho_np`, `epsilon_test` call functions that
      no longer exist in `bsk` (`epsilon`, `epsilon_tau`, `epsilon_delta_rho`, `g_e_*`) and
      raise `AttributeError`. Fix, make private, or delete.
- [ ] **Docstring checklist** for every public function: one-line summary · formula · Args with
      units and "per species, both spin components" where relevant · "accepts scalars and
      arrays" · Returns with units · `:cite:` reference · "See also" with working links.
      Leftovers to clean: `epsilon_np` ("what is kappa? (no Eq.9 in Ref.41)", "TO DO"),
      `vcritical` (formula labelled v_L), `vLandau` (k_F "[MeV]"), `Meff_hydro` (see TODO.md).
- [ ] **Executable examples:** enable `sphinx.ext.doctest`, put `>>>` examples in key
      docstrings and in the Quickstart, run `make doctest` in CI.
- [ ] Run `make coverage` (extension already loaded) to list undocumented public functions.

---

## D4 — New content

- [ ] **Pairing models page:** the three asymmetric-matter interpolation schemes (Chamel 2009,
      BSkG3 Eq. 3, BSkG4 Eq. 6) with formulas, references, and a generated figure comparing
      Δ_n, Δ_p at several asymmetries. Decide on "Pairing P" (fill or remove).
- [ ] **Functional comparison** (after `TODO_edf.md` E7): parameter table **generated from the
      parametrization registry** (a script in `docs/source/` writing an `.rst` include, same
      rule as figures), plus E/A (SNM & NeuM), M*/M and Δ_n curves for all functionals.
- [ ] **Inner-crust tables:** move the data to CSV files (`docs/data/*.csv`, loaded with
      `csv-table :file:`), cite the source for each table, and state the functional used.
- [ ] **Figure gallery** page: every generated figure with a link to its script.
- [ ] **Changelog** (`CHANGELOG.md`, included in the docs) — there is none yet; start it at the
      0.2.0 release (`TODO_edf.md` E8).
- [ ] **How to cite:** BibTeX for libnest (from `CITATION.cff`) plus the papers to cite for
      each functional. Also fix README "References": it cites PRC 88 024308 (BSk22–26) as
      "the BSk functional", but the implemented one is BSk31 (PRC 93 034337).
- [ ] Pasta page: keep hidden until `libnest.pasta` has real functions (today `volume()`
      returns 0).

---

## D5 — Polish (any time, optional)

- [ ] Theme: `sphinx_rtd_theme` is fine; `furo` or `pydata-sphinx-theme` would add dark mode
      and better mobile layout.
- [ ] "Edit on GitHub" links (`html_context` with github user/repo/branch), `html_title`,
      favicon/logo.
- [ ] Versioned docs (stable/latest) — only once there are several release lines users rely on.
- [ ] Keep README and docs landing page in sync (same Quickstart snippet, ideally included from
      one file).

---

## Suggested order

| Step | What | Effort | Depends on |
|------|------|--------|------------|
| 1 | D0 page/figure fixes (inner-crust table, Landau figure, TODO boxes, topsecret, bib keys, installation page, typos) | ~1 day | — |
| 2 | D0 physics-statement fixes (`rho2tau`, `rho2kf` convention, `units` tables, M*/M units) | ~0.5 day | your sign-off |
| 3 | D1 infrastructure, then turn on `-W -n` in CI | 0.5–1 day | step 1 |
| 4 | D2 restructure (landing page, conventions page, user guide skeleton) | 1–2 days | step 3 |
| 5 | D3 API reference pass | 1–2 days | `TODO_edf.md` E2 |
| 6 | D4 content (pairing page, functional comparison) | 1–2 days | `TODO_edf.md` E4–E7 |
| 7 | D5 polish | as wanted | — |

## Decisions (2026-10-02)
1. ✅ Physics fixes confirmed: `rho2tau` with (3π²)^(2/3); `rho2kf` = one isospin, both spins.
2. ✅ Remove hard-coded paths (done). `real_data_plots` stays in the docs for now (D2 decides
   where).
3. Keep the Makefile figure pipeline for now; the `.. plot::` directive can be revisited later.
4. Keep `sphinx_rtd_theme` for now.
