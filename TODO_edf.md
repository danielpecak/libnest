# TODO_edf — Multiple energy-density-functional parametrizations

Goal: libnest computes uniform-matter quantities (E/A, pressure, effective masses, mean
fields, pairing) for **any** Skyrme-type parametrization — BSk16/22/24/25/31, BSkG3, BSkG4,
SLy4, SkM*, … — with one implementation of the formulas, an explicit choice of functional
at the call site, and every parametrization validated against published numbers.

This file replaces/expands `TODO.md` → P4. Companion plan: [`TODO_doc.md`](TODO_doc.md).

---

## 1. Current state — what blocks it

- **Parameters are module globals** in `libnest/bsk.py` (`T0..T5`, `X0..X5`, `T2X2`,
  `ALPHA/BETA/GAMMA`, …) and all 44 functions read them directly.
- **`libnest/bskg.py` is a placeholder for a contributor** to fill with other functionals
  (the BSkG family). Today it is a byte-identical copy of `bsk.py` (commit `44c1807`), so
  filling it by editing the copied formulas would fork the implementation — §2a shows how it
  fits the new design instead.
- **Functional-specific values hard-coded outside the parameter block:**
  - kinetic term of `energy_per_nucleon` uses `HBARC**2/MN`, `HBARC**2/MP` (bare masses),
    whereas BSkG3/BSkG4/SLy4/SkM* use one common ħ²/2m (20.7355…, see `bsk_constants.h`);
  - pairing cutoff ε_Λ = **6.5 MeV hard-coded** in `I()` (`bsk.py` ~1470); the C header lists
    E_cut = 7.961 MeV (BSkG3) and 7.919 MeV (BSkG4);
  - reference-gap fits Δ_NeuM(k_F), Δ_SM(k_F) (coefficients 3.37968…, 11.5586…, cutoffs
    k_F > 1.38 / 1.31) hard-coded in `neutron_pairing_field` / `symmetric_pairing_field`;
  - asymmetric-matter interpolation: `I()`/`v_pi` always use `neutron_ref_pairing_field`
    (Chamel 2009); `*_eq2` is a duplicate of it; `*_eq3` (BSkG3) and `*_eq6` (BSkG4) exist
    but cannot be selected.
- **Defined but unused:** `FNP, FNM, FPP, FPM, KAPPAN, KAPPAP, YW` — no function reads them
  (`epsilon_pi_np` takes κ as an argument).
- **Import cycle `definitions` ↔ `bsk`:** `definitions.mu_q`, `xiBCS`, `E_minigap_rho_n`
  import `bsk` lazily, and `bsk` imports `mu_q` from `definitions`. `mu_q` depends on the
  effective mass, i.e. on the functional, so it belongs to the functional layer.
- **Name clash:** `units.ALPHA` (fine-structure constant) vs `bsk.ALPHA` (Skyrme exponent).
- **Consumers to migrate:** `plots.py` (~20 calls), `real_data_plots.py` (~10),
  `definitions.py` (3), docs figure scripts (8), `tests/`, `main.py`, `examples/`.

Reference design already exists: `~/projects/bsk/hpc-engine/bsk/bsk_constants.h` holds
**10 parametrizations** (BSk16, BSk22, BSk24, BSk25, BSk31, BSkG3, BSkG4, SLy4, SkM*, t0t3)
behind **one** `bsk_functional.c`, with a per-force pairing-interpolation switch. Mirror that.

---

## 2. Decision: separate modules, or something else?

**Recommendation: one implementation + parametrizations as data + a functional object.**
A parametrization is *data*, not code: BSk16…BSkG4…SLy4 differ only in numbers and a few
switches (pairing scheme, cutoff, ħ²/2m convention). Separate modules are right for
*code organization* (`skyrme.py`, `pairing.py`, `parametrizations.py`), wrong for
*parametrizations*.

| Option | For | Against |
|---|---|---|
| **A. One module per parametrization** (`bsk31.py`, `bsk24.py`, `bskg3.py`, … each with its own copy of the formulas) | Trivial to start | 10 × ~1000 lines of copies; every bug fix N times (see §5 — there are several); copies drift; can't loop over functionals |
| **B. Global switch** `set_functional("BSk24")` (old `TODO.md` P4 suggestion) | Minimal diff; old calls keep working | Hidden global state: comparing two functionals in one script means toggling; results depend on call order; tests can leak state into each other; not thread-safe |
| **C. `functional=` keyword on every function** | Explicit | ~50 signatures; must be threaded through internal chains (`v_pi → I → mu_q → effMn → B_q`) — forget one and BSk31 silently leaks into a BSk24 result |
| **D. ✅ Functional object** — immutable parameters + one class holding the formulas | Explicit; internal calls go through `self`, so functionals cannot mix; comparing is a loop; user-defined sets are easy; old `bsk` API kept as a BSk31 instance | One-time refactor (mechanical: `T0` → `p.t0`) |

What it looks like for a user:

```python
from libnest.edf import get_functional, available_functionals

for name in ["BSk24", "BSk31", "BSkG4"]:
    f = get_functional(name)                 # case-insensitive: "bsk24" works
    print(name, f.energy_per_nucleon(0.08, 0.08), f.effMn(0.08, 0.0))

f = get_functional("BSkG4")
f.neutron_ref_pairing_field(0.05, 0.01)       # uses BSkG4's own scheme (Eq. 6)
f.with_pairing(scheme="chamel2009")           # new object, same Skyrme part
f.params.t0, f.params.reference

import libnest.bsk as bsk                     # unchanged: BSk31, module-level functions
bsk.energy_per_nucleon(0.08, 0.08); bsk.T0
```

### 2a. Family modules: `bsk.py` and `bskg.py` (decided 2026-10-02)

Both stay as the user-facing entry points of their family; **neither holds formulas** — those
live once in `edf/skyrme.py`. A family module holds the family's parameter sets (data) and a
module-level API bound to a default functional.

- `bsk.py` — BSk family: module-level functions = BSk31 (today's API, unchanged); parameter
  sets `BSK16, BSK22, BSK24, BSK25, BSK31`.
- `bskg.py` — BSkG family, **the contributor's module**: parameter sets `BSKG3, BSKG4`
  (Skyrme part + pairing part: Eq. 3 / Eq. 6 scheme, common ħ²/2m, E_cut, κ convention) with
  references; module-level functions bound to a default BSkG functional, like `bsk.py`.

```python
# libnest/bskg.py after E2 — adding a functional is one parameter block, not 1600 lines
from libnest.edf import SkyrmeParameters, PairingParameters, register

BSKG3 = register(
    SkyrmeParameters(name="BSkG3", t0=..., t1=..., ..., hbar2m_n=20.73553, hbar2m_p=20.73553,
                     reference="Grams et al., EPJA 59, 270 (2023)"),
    PairingParameters(scheme="bskg3_eq3", cutoff=7.961, ...))
BSKG4 = register(..., PairingParameters(scheme="bskg4_eq6", cutoff=7.919, ...))

_default = BSKG4                              # see open question 7
energy_per_nucleon = _default.energy_per_nucleon
...
```

`get_functional("BSkG4")` and `libnest.bskg.BSKG4` are the same object (the registry imports
the family modules lazily on first lookup).

---

## 3. Target design

```
libnest/
  edf/
    __init__.py          get_functional, available_functionals, register, SkyrmeFunctional
    parameters.py        SkyrmeParameters, PairingParameters (frozen dataclasses)
    parametrizations.py  registry (get_functional, register); standard Skyrme sets SLy4,
                         SkM*, t0t3 — family sets live in their family modules
    skyrme.py            SkyrmeFunctional: E/A, pressure, c_s, B_q, U_q, M*, C^ρ, C^τ, ε_*, v_π, I
    pairing.py           gap models (Δ_NeuM, Δ_SM fits) + interpolation schemes registry
  bsk.py                 BSk family: BSk16/22/24/25/31 sets; module-level API = BSk31
                         (compatibility: old names & constants)
  bskg.py                BSkG family (contributor): BSkG3/BSkG4 sets; module-level API
  definitions.py         functional-independent only (rho2kf, eF_n, vLandau, …)
```

```python
@dataclass(frozen=True)
class SkyrmeParameters:
    name: str
    t0: float; t1: float; t2: float; t3: float
    x0: float; x1: float; t2x2: float; x3: float      # t2*x2 stored directly (BSk: t2 = 0)
    alpha: float
    t4: float = 0.0; t5: float = 0.0; x4: float = 0.0; x5: float = 0.0
    beta: float = 0.0; gamma: float = 0.0              # standard Skyrme: t4 = t5 = 0
    hbar2m_n: float = HBAR2M_n; hbar2m_p: float = HBAR2M_p
    # metadata (not used in uniform matter): spin-orbit, Wigner, odd/even pairing factors
    w0: Optional[float] = None; w0_prime: Optional[float] = None
    reference: str = ""; doi: str = ""; notes: str = ""

@dataclass(frozen=True)
class PairingParameters:
    scheme: str = "chamel2009"        # | "bskg3_eq3" | "bskg4_eq6"
    gap_model: str = "cao2006"        # which Δ_NeuM / Δ_SM fit (today's hard-coded one)
    cutoff: float = 6.5               # ε_Λ [MeV]
    kappa_n: float = 0.0; kappa_p: float = 0.0
    kappa_convention: str = "absolute"   # BSkG uses "relative" (see §6)
    f_np: float = 1.0; f_nm: float = 1.0; f_pp: float = 1.0; f_pm: float = 1.0
```

Rules:
- **Move first, rename later.** Methods keep today's names (`effMn`, `B_q`, `U_q`, …) so the
  refactor is mechanical and reviewable; nicer names can come later, with aliases.
- Keep the established numerics: `np.asarray` in, scalar out for scalar input, `DENSEPSILON`
  before every division by ρ, the scalar early-return pattern of the pairing fields.
- Functional-independent pieces stay plain functions: `Lambda(x)`, generic numerical
  derivatives (take a callable), everything in `definitions`.
- `definitions.mu_q / xiBCS / E_minigap_rho_n` become methods; `definitions` keeps thin
  wrappers with an optional `functional=` argument (default BSk31) → the import cycle is gone
  (`edf` imports `definitions`, never the reverse).
- `plots.py` functions get `functional="BSk31"` (name or object) — and comparison plots that
  take a list of functionals.

---

## 4. Parametrizations to port

Source of the numbers: `bsk_constants.h` (each block cites its paper; re-check every value
against the paper table when porting). "bask24" in the request = **BSk24**.

| Name | Family | Reference | Structural notes (from `bsk_constants.h`) | Pairing scheme | ħ²/2m |
|---|---|---|---|---|---|
| BSk16 | BSk | Chamel+ NPA 812 72 (2008) | t4 = t5 = 0, β = γ = 0, t2 ≠ 0, α = 0.3 | Chamel 2009 | per species |
| BSk22 | BSk | Goriely+ PRC 88 024308 (2013) | α = γ = 1/12, β = 1/2; t2x2 recovered from MOCCa decomposition; κ = 0 | Chamel 2009 | per species |
| **BSk24** | BSk | Goriely+ PRC 88 024308 (2013) | same as BSk22 | Chamel 2009 | per species |
| BSk25 | BSk | Goriely+ PRC 88 024308 (2013) | same as BSk22 | Chamel 2009 | per species |
| BSk31 | BSk | Goriely+ PRC 93 034337 (2016) | existing; α = 1/5, β = 1/12, γ = 1/4; κ_n, κ_p absolute | Chamel 2009 | per species |
| **BSkG3** | BSkG | Grams+ EPJA 59 270 (2023) | exponents as BSk31; t2 = 0.01; κ relative (disabled in C); E_cut 7.961 MeV | **Eq. 3** | common 20.73553 |
| **BSkG4** | BSkG | Grams+ arXiv:2411.08007 (2024) | as BSkG3; E_cut 7.919 MeV | **Eq. 6** | common 20.73553 |
| SLy4 | Skyrme | Chabanat+ NPA 635 231 (1998) | standard Skyrme: t4 = t5 = 0, α = 1/6 | Chamel 2009 (gap model only) | common |
| SkM* | Skyrme | Bartel+ NPA 386 79 (1982) | standard Skyrme, α = 1/6 | Chamel 2009 (gap model only) | common |
| t0t3 | toy | Saclay validation package | only t0, t3; α = 1 → closed-form E/A, ideal test case | — | common |

Note: the C BSk16 block carries BSk31's κ and f^± values (looks copy-pasted) — check
against Chamel 2008 before porting.

---

## 5. Physics bugs found during this analysis — fix once, in the new class

These must be fixed **before** porting new parametrizations (phase E3), otherwise the new
functionals get validated against buggy code. Each fix: own commit, invariant test, your
sign-off. **Status: all done in E3 (2026-10-02)** — see the commits listed under E3.

- [x] **`mu_q` treated M*/M as a mass in MeV** (31096 MeV at ρ_n = 0.08) → now
      μ_q = ħ²k_F²/2M*_q = B_q k_F² = 33.1 MeV. Also fixed: `mu_q(array, 0.0, 'p')` raised
      `ValueError` (`np.where` on a 0-d array, NumPy 2).
- [x] **`v_pi` had the same unit error** (−1.3e5 MeV fm³) → B_q^{3/2}; BSk31 neutron matter
      now −714 … −213 MeV fm³ for ρ_n = 0.005 … 0.05 fm⁻³.
- [x] **`U_q` t3 term** → exact derivative of `epsilon_rho_np` (as hpc-engine `g_U_rho_n`
      computes; its doc comment has the same garbled formula, its code is right). Docstring
      states that τ- and gradient-dependent parts are not included.
- [x] **Gradient terms:** `epsilon_np` uses (∇ρ_n + ∇ρ_p)²; the cross term in
      `epsilon_delta_rho_np` is (|∇ρ|² − |∇ρ_n|² − |∇ρ_p|²)/2.
- [x] **`rho2tau` factor π^(2/3)** — fixed in `TODO_doc.md` D0.
- [x] **Kinetic term** of `energy_per_nucleon` and `epsilon_np` uses the parameter set's
      ħ²/2m (BSk31 changes ≤ 2e-7 relative).
- [x] **BSkG4 Eq. 6 floor — analysed, deliberately not changed.** Over a 641 × 201 grid
      (ρ ≤ 0.32 fm⁻³, all proton fractions) the C-style 1e-8 floor changes Δ_n, Δ_p by at most
      9e-7 MeV, and only where the gap is closed anyway (no point with a gap > 1 keV differs).
      The E4 parity tests against C therefore use `atol ≈ 1e-6 MeV` instead.
- [x] **`testMe`** removed; its closed-form symmetric-matter E/A (PRC 80 065804 Eq. A13, with
      the n/p-averaged ħ²/2m) is now a test of the general formula. **`plots.py`:** three
      energy-density plots rewired to the renamed `bsk` functions, four dead ones removed.

---

## 6. Migration plan (each phase = one PR, full test suite green after each)

### E0 — Safety net (no code changes yet)  ✅ DONE (branch `edf-refactor`, 2026-10-02)
- [x] **Keep `libnest/bskg.py` untouched** — it is the contributor's placeholder; E2 turns it
      into the BSkG family module (§2a).
- [x] **Golden-master test for BSk31.** `tests/golden/spec.py` (cases + input grid),
      `tests/golden/make_golden.py` (generator; without `--write` it only reports which
      cases would change), `tests/golden/bsk31.json` (inputs, results, API snapshot — JSON,
      one value per line, so a git diff shows every changed number),
      `tests/test_golden_bsk31.py`.
      - 51 cases: all 43 public `bsk` functions (`U_q`, `B_q`, `v_pi`, `I` for both q) plus
        `definitions.mu_q`, `xiBCS`, `E_minigap_rho_n`, which move in E2;
      - grid: 20 total densities × 5 proton fractions (0 … 0.7: NeuM, SNM, proton-rich),
        including ρ = 0 and points on both sides of every pairing cutoff; array input and
        point-by-point Python-float input;
      - current edge behaviour is recorded, so the refactor must keep it (or change it on
        purpose): scalar ρ = 0 raises `ZeroDivisionError` in the 6 analytic neutron-matter
        functions (array input gives `inf`); `I`, `v_pi`, `epsilon_pi_np`, `epsilon_np` return
        masked arrays;
      - tolerance `rtol=1e-12`, `atol=1e-15`;
      - mutation-checked: T5 changed by 1e-4 % → 34 cases fail; a pairing cutoff moved by
        0.005 fm⁻¹ (scalar or array branch) → 12 fail; Λ(x) constant changed by 1e-7 → 7 fail;
        `testMe` renamed → API test fails.
- [x] **Public API snapshot** in the same file: parameter lists of all public `bsk` functions
      and of the three `definitions` helpers, values of `T0 … KAPPAP`. The test accepts new
      parameters only when appended with defaults.

### E1 — Parameters become data (no behavior change)  ✅ DONE (2026-10-02)
- [x] `libnest/edf/parameters.py`: frozen dataclasses `SkyrmeParameters` (t0…t5, x0…x5,
      t2x2, α/β/γ, ħ²/2m per species, `yw`, reference/doi/notes) and `PairingParameters`
      (cutoff ε_Λ, κ_n/κ_p, f^±). Fields not used before E4 (`scheme`, `gap_model`,
      `kappa_convention`) are deliberately **not** added yet, so nothing can be set and
      silently ignored.
- [x] The BSk31 values live in `bsk.py` (`BSK31 = register(SkyrmeParameters(...),
      PairingParameters(...))`); the old constants `T0 … KAPPAP` are read from it.

### E2 — One implementation (no behavior change)  ✅ DONE (2026-10-02)
- [x] `libnest/edf/skyrme.py`: `SkyrmeFunctional` with all 42 functional-dependent
      functions as methods, plus `mu_q` moved from `definitions`. Generated mechanically
      from the old `bsk.py` (AST-based: `T0` → `p.t0` and internal calls → `self.` in code
      lines only, docstrings untouched); every transformed line was reviewed. `Lambda` and
      `_clamp_nonnegative` stay plain functions. `HBAR2M_n/p` now come from the parameter
      set (same values); the pairing cutoff 6.5 MeV in `I()` is `self.pairing.cutoff`.
      Bare masses (`HBARC**2/MN` in `energy_per_nucleon`, `epsilon_np`) unchanged → E3.
- [x] `libnest/edf/parametrizations.py`: `register`, `get_functional` (case-insensitive),
      `available_functionals`; family modules (`libnest.bsk`, `libnest.bskg`) are imported
      lazily on first lookup. `libnest/edf/__init__.py` re-exports the public API.
- [x] `bsk.py` = BSk family module: `BSK31`, the old constants, module-level functions =
      bound methods of `BSK31`, `__all__` (so the docs still list them as `bsk` functions),
      old re-exports (`bsk.MN`, `bsk.rho2kf`, `bsk.mu_q`, …) kept.
- [x] `definitions.mu_q`, `xiBCS`, `E_minigap_rho_n` take an optional `functional=` (object
      or name, default BSk31) and import `edf` lazily → no module-level import cycle.
      _Deviation from the plan:_ `xiBCS` and `E_minigap_rho_n` stay in `definitions` (they
      only need a pairing field); only `mu_q` became a method.
- [x] **Verified:** golden master bit-identical in all 51 cases and every field (array,
      scalar, masks, scalar exceptions) — `bsk31.json` untouched; API snapshot compatible
      (only the appended `functional=None`); 88 tests green (new `tests/test_edf.py`:
      registry, immutability, independence of functionals, `functional=` wrappers, cutoff
      wiring, scalar/array); `test_imports` now walks subpackages; ruff clean; strict docs
      build (`-W -n`) green with a new `docs/edf.rst` page.
- [x] The commented-out legacy NeuM block at the end of the old `bsk.py` (`C_rho`,
      `epsilon_rho`, … — dead code) was dropped; it remains in git history.
- `libnest/bskg.py` untouched (still the old copy); it imports fine and registers nothing.

### E3 — Fix the physics bugs from §5 (numbers change on purpose)  ✅ DONE (2026-10-02)
One commit per fix: invariant test first (checked to fail on the old code) → fix →
`make_golden` (lists the affected cases) → `--write`; the reason is in each commit message
and the git diff of `bsk31.json` shows exactly which numbers moved.

| Commit | Fix | Golden cases changed | Invariant test (`tests/test_skyrme.py`) |
|---|---|---|---|
| `aa35294` | `mu_q` units | mu_q, I, v_pi, ε_π, ε | free gas: μ = e_F; μ/e_F = M/M* |
| `68eaf58` | `mu_q` mixed scalar/array input | — | array ρ_n with scalar ρ_p, both q |
| `7749e9c` | `v_pi` = −8π²/I · B_q^{3/2} | v_pi, ε_π, ε | gap-equation normalization; attractive, < 2000 MeV fm³ |
| `f41feba` | `U_q` t3 term | U_q | U_q = ∂ε_ρ/∂ρ_q; B_q = ħ²/2m + ∂ε_τ/∂τ_q |
| `bee7e0e` | ε: (∇ρ_n + ∇ρ_p)² | ε | gradient part of ε = ε_Δρ with the total gradient |
| `f081c5e` | ε_Δρ cross term | ε_Δρ, ε | n ↔ p exchange symmetry |
| `12b41f5` | kinetic term from ħ²/2m | E/A and everything built on it, ε | free gas with custom ħ²/2m; general vs analytic NeuM E/A |
| `ae738b6` | remove `testMe` | testMe case + API entry | closed-form SNM limit of Eq. (A13) |
| `6c40b4a` | `plots.py` dead functions | — | `tests/test_plots.py`: every plot runs and renders |

Suite after E3: 104 tests + 181 subtests, ruff clean, strict docs build green.

### E4 — Pairing as a component
- [ ] `edf/pairing.py`: gap model (Δ_NeuM, Δ_SM fits + k_F cutoffs) and a scheme registry
      {`chamel2009` (= today's `_ref` and `_eq2`), `bskg3_eq3`, `bskg4_eq6`}; each functional
      has a default scheme, overridable with `with_pairing(scheme=...)`.
- [ ] ε_Λ from `PairingParameters.cutoff` instead of the hard-coded 6.5 MeV.
- [ ] Old function names (`neutron_ref_pairing_field_eq3`, …) stay in `bsk` as wrappers.
- [ ] Parity test against the C `g_Delta_n_ref` / `g_Delta_p_ref` for all three schemes.

### E5 — Port the parametrizations (§4)
- [ ] One registry entry per force with reference, DOI and the notes from `bsk_constants.h`.
- [ ] Case-insensitive lookup + aliases (`"bsk24"`, `"SkMstar"`); `register()` for
      user-defined sets (`dataclasses.replace(BSK31, t0=...)` for sensitivity studies).

### E6 — Validation for every functional
- [ ] **Universal invariants**, parametrized over the registry: E/A → 0 for ρ → 0; SNM bound
      at saturation; NeuM above SNM; M* > 0; scalar result == array result; no NaN at ρ = 0
      (watch `ρ^(β−1)` with β = 0 — fine only because of `DENSEPSILON`).
- [ ] **Published saturation properties** per force — n₀, a_v, J, L, K_v, M*_s/M — stored as
      `reference_values` in the registry entry with table/page cited; tests read them from
      there, so adding a functional is self-validating. (Don't type them from memory —
      copy from the paper tables.)
- [ ] **t0t3:** closed-form E/A check.
- [ ] **Cross-code parity with hpc-engine:** generate reference tables (E/A, B_q, Δ_n, Δ_p on
      a grid) with the C uniform-matter tools for every `BSK` id and compare. Note
      `HBARC` differs (libnest 197.3269804 vs C 197.32697881, ~1e-8 relative) — align or set
      the tolerance accordingly.

### E7 — Consumers and docs
- [ ] `plots.py`: `functional=` argument + comparison plots; `real_data_plots.py`: same
      argument, default BSk31.
- [ ] `main.py`, README Quick Start, `examples/`: show `get_functional`.
- [ ] Docs: parameter table **generated from the registry** (script in `docs/source/`, same
      rule as figures — replaces the hand-written BSk31 table in the `bsk` docstring);
      comparison figures; a "Choosing a functional" page (see `TODO_doc.md` D4).

### E8 — Release
- [ ] `CHANGELOG.md`, version 0.2.0 (new API added, old API intact).
- [ ] Point `TODO.md` P4 to this file; update `CLAUDE.md` architecture (`edf` layer between
      `definitions` and `bsk`).

---

## 6a. Contributor track — `bskg.py`

Work that does not depend on the refactor and can start now (it feeds E5/E6):
- [ ] BSkG3 / BSkG4 parameter tables from the papers (Grams+ EPJA 59 270 (2023);
      Grams+ arXiv:2411.08007), cross-checked against `bsk_constants.h`, including ħ²/2m,
      W0, W0′.
- [ ] Published infinite-matter properties of both forces (n₀, a_v, J, L, K_v, M*_s/M,
      M*_v/M) with table/page — they become the E6 reference values.
- [ ] Pairing: confirm the Eq. 3 (BSkG3) / Eq. 6 (BSkG4) schemes already implemented in
      `bsk.py` (`*_eq3`, `*_eq6`), the cutoff E_cut and the relative-κ convention (question 2).
- [ ] Reference curves from the papers (E/A in SNM and NeuM, gaps) for later comparison
      figures.

Please **don't edit the copied formulas in `bskg.py`** before E2 — they will be replaced;
formula fixes go to `bsk.py` (the one implementation). After E2, `bskg.py` is filled as in
§2a: one parameter block per functional.

---

## 7. Open questions (physics decisions — yours)

1. **Reference gaps per functional?** The C code uses the same Δ_NeuM/Δ_SM fit for every
   force. `physics.rst` says BSk16–17 were fitted to gaps *without* self-energy and BSk30–32
   *with* it. Which curve is the current fit (NeST.pdf Eqs. 5.11–5.12), and should BSk22–25
   use a different one? Proposal: start with the C behavior (parity) and keep `gap_model`
   as the switch for later.
2. **BSkG relative κ** (κ_n = 123.20, κ_p = 129.07 fm⁸ with g_q = V_q[1 + κ_q(∇ρ)²]):
   implement that convention properly, or keep it disabled as in C?
3. **f^± factors:** if, as I understand, they distinguish even/odd nucleon numbers in finite
   nuclei, they're metadata in uniform matter (f⁺ = 1). Confirm, or tell me where they enter.
4. ✅ **Package name:** `libnest.edf` (default applied in E1–E2, 2026-10-02).
5. ✅ **Future of `libnest.bsk`:** stays permanently as the BSk family module / BSk31 shortcut
   (default applied in E2, 2026-10-02).
6. **Single source of truth with hpc-engine:** later, both `bsk_constants.h` and
   `parametrizations.py` could be generated from one data file (e.g. TOML). Worth it, or are
   parity tests enough?
7. **Default of `libnest.bskg`'s module-level functions:** BSkG4 (newest, proposed) or BSkG3?

## 8. Effort estimate

| Phase | Effort |
|---|---|
| E0 safety net | 0.5 day |
| E1 parameters as data | 0.5 day |
| E2 one implementation + compat layer | 1–2 days |
| E3 physics bug fixes | 1–2 days (+ checking derivations) |
| E4 pairing component | 1 day |
| E5 port 9 parametrizations | 1 day |
| E6 validation | 1–2 days |
| E7 consumers + docs | 1–2 days |
| E8 release | 0.5 day |

Out of scope: finite nuclei, spin-orbit, Wigner and Coulomb terms (uniform matter only —
W0/W0' are stored as metadata).
