# TODO — Roadmap to a pip-installable, student-friendly libnest

Goal: a student can run `pip install libnest`, then in a few lines import it, analyze
their simulation data, and benchmark it against known physics — without touching the repo
internals or editing hard-coded paths.

Priorities: **P0** = blocks *everything* (library is unusable until done) · **P1** = blocks
`pip install` · **P2** = blocks "analyze / benchmark my data" · **P3** = polish & trust.

---

## P0 — Make the library import and pass its tests  ✅ DONE

Nothing else mattered until `import libnest.bsk` worked. It now does.

- [x] **Fixed the `IndentationError`s in `libnest/bsk.py`.** Three over-indented bodies:
      `neutron_ref_pairing_field` (~305), `proton_ref_pairing_field` (~340), and
      `epsilon_delta_rho_np` (~1001). All compile and import now.
- [x] **Fixed a deeper scalar-path bug surfaced while testing:** `neutron_pairing_field` /
      `symmetric_pairing_field` called `np.where()` on a 0-d array before the scalar
      early-return, which NumPy 2.x rejects — every scalar call crashed. Moved the `np.where`
      into the array-only branch.
- [x] **Added `tests/test_imports.py`** — imports every `libnest` submodule (auto-discovered).
- [x] **Added regression tests** in `tests/test_bsk.py` for the reference pairing fields, the
      gradient term, and scalar input. _Full suite: 34 tests + 13 subtests green in a clean venv._
- [x] **Made `main.py` a CI-checked smoke test** — the `tests` workflow now runs
      `python main.py` end-to-end on every push, alongside the import + regression tests.

## P1 — Make it installable from PyPI  ✅ DONE (except the actual publish)

A student's entry point is `pip install`. The packaging now exists and is validated.

- [x] **Added `pyproject.toml`** (PEP 621, setuptools backend). Version single-sourced from
      `libnest.__version__`; Python `>=3.9`; runtime deps `numpy/scipy/matplotlib/pandas`;
      `test` and `docs` extras; project URLs; pytest config.
      _Verified:_ `pip install -e ".[test]"` works in a clean venv; `python -m build` +
      `twine check` both PASS for the sdist and wheel.
- [x] **Added `LICENSE`** (MIT, `2022-2026 Daniel Pęcak`). ⚠️ **Confirm the author is happy
      with MIT and the copyright line** — it matches the README's existing claim but should be
      signed off.
- [x] **Split dependencies.** `requirements.txt` is now runtime-only; docs deps live in
      `docs/requirements.txt` / the `[docs]` extra; tests in the `[test]` extra.
- [x] **Pinned floors** (`numpy>=1.23`, `pandas>=2.0`, `scipy>=1.9`, `matplotlib>=3.6`) to
      avoid the NumPy-2/pandas ABI import failure seen on the system Python.
- [x] **Added `.github/workflows/release.yml`** — builds sdist+wheel on `v*` tags, runs
      `twine check`, and publishes to PyPI via Trusted Publishing (OIDC).
- [x] **Published.** Dry-run to TestPyPI first (via a manual `workflow_dispatch` job), then
      tagged `v0.1.0` → **`libnest 0.1.0` is live on PyPI** (`pip install libnest` works).
- [x] **Also added `.github/workflows/tests.yml`** (pytest + coverage on Python 3.9–3.12) —
      this was a P3 item; done early so packaging changes stay green.

## P2 — Make it usable for *my own data* (analyze + benchmark)  ✅ DONE

The scientific workflow a student actually wants. The data path is now portable, the reader
works on modern pandas, and there is a runnable analyze-your-data example + physics benchmarks.

- [x] **Fixed `libnest/myio.py` (data loader).** `txt2df` now uses `df = df.set_axis(...)`
      (pandas ≥ 2.0) and a raw docstring; `readDimTxt` no longer calls an unimported `sys.exit`,
      defaults missing dimensions to 1, and works for 1D/2D/3D (removed the hard 2D exit).
      Verified by `tests/test_myio.py` (9 tests, incl. 3D + missing-dim regressions).
- [x] **Removed hard-coded absolute paths from `libnest/real_data_plots.py`.** Now driven by
      `DATA_ROOT = os.environ.get("LIBNEST_DATA", ".")` — portable, overridable per user.
      (Follow-up, still open: pass the data dir as a function *argument* rather than a module
      global — see P3.)
- [x] **Example dataset is generated at runtime, not stored.** `examples/analyze_data.py`
      writes a tiny 8×8 density map on the fly (`generate_sample_dataset`), so no data files
      live in git. `.gitignore` keeps generated `*/data/*.txt` and the output PNG out of the repo.
- [x] **Wrote `examples/analyze_data.py`** — the copy-paste starting point: generate/load a
      density map via `myio`, compute per-cell `kF` and neutron pairing gap via
      `definitions`/`bsk`, print a summary, and save a 2D gap-map figure. Point `LIBNEST_DATA`
      at your own data (same format) to run it for real.
- [x] **Added `tests/test_benchmark.py`** (5 tests) asserting against published BSk31 values
      with cited sources: saturation density n0 ≈ 0.1587 fm⁻³, E/A ≈ −16.05 MeV, SNM/neutron
      Fermi momenta (1.33 / 1.68 fm⁻¹), and the O(1 MeV) neutron pairing-gap scale. Doubles as
      regression protection for the functional.

## P3 — Coherence, trust, and polish

- [x] **Unify the test suite.** Deleted the print-script `tests/tests.py` and the
      never-discovered `tests/units_tests.py`; their intent lives in `tests/test_units.py`
      (now pytest-discoverable). Deleted the dead `examples/legacy_tests.py` scratchpad and
      pointed the README at `examples/analyze_data.py` instead. Added `tests/test_tools.py`
      (8 tests) — which surfaced and fixed a crash in `threeSlice` on 2D input (`np.asarray`
      was given two positional args). Suite now: 64 tests + 17 subtests.
      _Still open:_ wire a coverage badge (`pytest-cov` is already in the `[test]` extra / CI).
- [x] **Add a CI test workflow** (`.github/workflows/tests.yml`) running pytest + coverage on
      3.9–3.12. _Done as part of P1._
- [ ] **Resolve or file the ~21 in-code TODOs** catalogued in `docs/history/ANALIZA.md` (missing imports,
      duplicated `E_minigap_*` functions, `nucleus` background-density handling). Each should
      become a fixed bug or a tracked issue, not a silent comment.
- [x] **Audited cross-module dependencies** with `pyflakes`. `tools.py` imports and
      `definitions.mu_q` are clean (no undefined names). The audit turned up a *real* latent
      bug: `definitions.Meff_hydro` had a copy-pasted body `return HBAR2M_n * kF**2` using an
      undefined `kF` — replaced with the documented formula and covered by tests. Also fixed 4
      invalid-escape `SyntaxWarning`s (raw strings in `plots`/`tools`/`units` docstrings).
      _Still open (cosmetic):_ ~11 unused imports and 3 unused locals (`ax`, `V`).
      ⚠️ **Physics to confirm:** `Meff_hydro`'s docstring says "in units neutron mass" but its
      LaTeX includes `m_n`; the two are inconsistent — verify against Magierski (2004).
      ⚠️ **Physics to confirm:** `tools.condensationEnergy` ignores its documented `dV`
      argument, and the code divides by `rho` (`|Δ|²/(ε*_F·ρ)`) while the docstring formula
      multiplies by it (`|Δ|²/ε*_F·ρ`). Left as-is (only smoke-tested) pending your check.
- [~] **Linter/formatter + type hints.** `ruff` is wired as a **linter** (`select = ["F"]`,
      `docs/` excluded, `F841` ignored), enforced by a `lint` job in CI; a `[lint]` extra
      installs it. This removed 12 unused imports. _Deliberately not done:_ `ruff format`
      (wholesale reformat — big diff, skipped by choice). _Still open:_ type hints and a
      pre-commit hook.
- [x] **README cleanup.** Added PyPI-version + Tests-workflow badges; replaced the dead
      ReadTheDocs link with the GitHub Pages docs URL (`danielpecak.github.io/libnest`).
      (Python badge 3.9+, version/"Last Updated", venv/`pip install` section already done.)
      _Note:_ the Pages link needs GitHub → Settings → Pages set to serve from the `sphinx`
      branch (where `documentation.yml` deploys) — otherwise it 404s.
- [x] **Consolidated the Polish planning docs.** Moved `ANALIZA.md` and
      `PODSUMOWANIE_ZMIAN.md` to `docs/history/` with a "historical snapshot" header; `TODO.md`
      is the single living roadmap.

## P4 — Support Skyrme functionals beyond BSk31

Right now `libnest/bsk.py` hard-codes the BSk31 parameters as module-level globals
(`T0..T5`, `X0..X5`, `T2X2`, `ALPHA/BETA/GAMMA`, `YW`, `FNP/FNM/FPP/FPM`,
`KAPPAN/KAPPAP`, `HBAR2M_n/HBAR2M_p`) and every function reads those globals. The C
reference `bsk_constants.h` (project `bsk/hpc-engine`) already organizes six
parametrizations behind one interface — mirror that design in Python.

- [ ] **Refactor to a selectable parameter set.** Replace the flat globals with one
      parameter set per functional (a `dataclass` or dict keyed by name, e.g.
      `FUNCTIONALS["BSk31"]`), and a way to pick the active one. Decide the API:
      - a module-level `set_functional("BSk31")` that swaps the active set (simplest,
        keeps existing call sites working), **or**
      - an optional `functional=` argument threaded through the public functions (purer,
        bigger diff). _Recommendation: `set_functional()` + a default of BSk31 for
        backward compatibility, so existing code and tests keep passing._
- [ ] **Port the six parametrizations** from `bsk_constants.h` with their references:
      `BSk16` (Chamel+ NPA 812, 2008), `BSk31` (existing; PRC 93, 2016),
      `BSkG3` (Grams+ EPJA 59, 2023), `BSkG4` (Grams+ arXiv:2411.08007, 2024),
      `SLy4` (Chabanat+ NPA 635, 1998), `SkM*` (Bartel+ NPA 386, 1982).
- [ ] **Honor the structural differences** (documented in `bsk_constants.h`):
      - standard Skyrme forces (`SLy4`, `SkM*`) have `T4=T5=0`, `BETA=GAMMA=0` — the
        extended-Skyrme gradient terms vanish; verify libnest's formulas reduce correctly;
      - `BSkG3/BSkG4` use a single `hbar^2/2m` for both species and set `KAPPAN=KAPPAP=0`
        (their gradient-pairing uses a different, *relative* convention — disabled here);
      - pairing effective-mass factors `FNP/FNM/FPP/FPM` differ per force (BSk uses ≠1,
        the others 1.0).
- [ ] **Pairing interpolation options.** Add the two asymmetric-matter schemes from the C
      code: `PAIRING_CHAMEL` (linear, PRL 102 152503 2009 — current default) and
      `PAIRING_BSKG4_EQ6` (geometric, arXiv:2411.08007 Eq. 6, default for BSkG4).
      `neutron_ref_pairing_field`/`proton_ref_pairing_field` already implement the Chamel
      ansatz; add the BSkG4 variant and select per functional.
- [ ] **Benchmark each functional.** Extend `tests/test_benchmark.py` with the published
      saturation point (n0, E/A) per force, so adding a parametrization is self-validating.
- [ ] **Docs.** Turn the hard-coded "BSk31 parameters" table in `bsk.py`'s docstring into a
      table per functional, and note how to select one (`set_functional`).

---

### Suggested order of attack
1. **P0** (half a day) → library works again, tests green.
2. **P1 packaging + LICENSE** (1 day) → `pip install libnest` works.
3. **P2 myio + real_data_plots + example data** (1–2 days) → students can load and plot
   their own data.
4. **P2 benchmark suite + P3 CI** (1–2 days) → results are trustworthy and stay that way.
5. **P3 polish** ongoing.
