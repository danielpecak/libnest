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
      global — tracked in P3.)
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
- [ ] **Resolve or file the in-code TODOs.** The ~21 of `docs/history/ANALIZA.md` are down to
      (2026-10-17): `nucleus.py` ×5 (background density for both species, 3D directions,
      formula/normalization of the documented quantity), `edf/skyrme.py` ×3 (two wish-list
      comments — M* plots, non-uniform matter — and the `epsilon_np` docstring "TO DO: write eq
      for jsum2 and jdiff2"), one hidden `.. todo::` in `docs/inner-crust.rst` (Landau
      velocity). `bskg.py` carries copies; it is the contributor's (`TODO_edf.md` §6a). Each
      should become a fixed bug or a tracked issue, not a silent comment.
- [ ] **`real_data_plots`: data directory as a function argument** instead of the module
      globals `TXT_PATH*` (follow-up of P2).
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

## P4 — Support Skyrme functionals beyond BSk31  → see [`TODO_edf.md`](TODO_edf.md)

Planned and tracked in `TODO_edf.md`. The design differs from the first idea here (a global
`set_functional()`): one implementation of the formulas, parameter sets as data, and a
functional object (`libnest.edf.get_functional("BSk31")`); `libnest.bsk` keeps its API as
the BSk31 instance.

Status (2026-10-17, branch `edf-refactor`): E0–E4 done — golden master, parameters as data,
`SkyrmeFunctional`, physics fixes, pairing schemes (Chamel 2009, BSkG3 Eq. 3, BSkG4 Eq. 6)
with C parity tests, relative κ of BSkG. Open: E5 port the parametrizations (BSkG3/BSkG4 by
the contributor in `bskg.py`), E6 validation, E7 consumers and docs, E8 release 0.2.0.

---

### Suggested order of attack
P0–P2 are done (`libnest 0.1.0` on PyPI). Current work: P4 = `TODO_edf.md` (E5–E8, release
0.2.0) and the documentation plan `TODO_doc.md` (D2–D4); P3 polish alongside.
