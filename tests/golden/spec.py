# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Daniel Pęcak
"""
Golden-master specification for the BSk31 functional (TODO_edf.md, phase E0).

Every case calls one public function on stored inputs. ``make_golden.py`` evaluates
the cases and writes ``bsk31.json``; ``tests/test_golden_bsk31.py`` evaluates them
again and compares. The refactor phases E1-E2 must not change any of these numbers;
intentional changes (E3 physics fixes) regenerate the file, so the git diff of
``bsk31.json`` shows exactly which values moved.
"""
import importlib
import inspect

import numpy as np

GOLDEN_FILE = "bsk31.json"


def make_inputs():
    """Returns the input arrays (dict of 1D numpy arrays) shared by all cases.

    The grid spans total density x proton fraction, including rho = 0, pure
    neutron matter, symmetric matter and proton-rich matter (eta < 0). The values
    are stored in the golden file, so the test does not depend on this function.
    """
    # Includes points on both sides of the pairing cutoffs: kF(rho_n) = 1.38 at
    # rho_n ~ 0.0888 (Delta_NeuM), kF(rho/2) = 1.31 at rho ~ 0.1519 (Delta_SM), and
    # kF(rho/2) = 1.38 at rho ~ 0.1775 (Delta_NeuM(rho/2) in the Eq. 6 scheme).
    rho_grid = [0., 1e-4, 1e-3, 0.01, 0.02, 0.04, 0.06, 0.08, 0.0885, 0.089, 0.0895,
                0.1, 0.12, 0.151, 0.152, 0.16, 0.177, 0.178, 0.2, 0.3]
    yp_grid = [0., 0.1, 0.3, 0.5, 0.7]
    rho, yp = (a.ravel() for a in np.meshgrid(rho_grid, yp_grid, indexing="ij"))
    rho_n = (1. - yp)*rho
    rho_p = yp*rho
    n = rho.size
    rho_grad_n = 0.5*rho_n          # [fm^-4], arbitrary but non-trivial
    rho_grad_p = -0.3*rho_p
    return {
        "rho": rho,
        "rho_n": rho_n,
        "rho_p": rho_p,
        "kF_n": (3.*np.pi**2*rho_n)**(1./3.),
        "tau_n": 0.6*(3.*np.pi**2*rho_n)**(2./3.)*rho_n,
        "tau_p": 0.6*(3.*np.pi**2*rho_p)**(2./3.)*rho_p,
        "jsum2": (0.1*rho)**2,
        "jdiff2": (0.05*(rho_n - rho_p))**2,
        "rho_grad_n": rho_grad_n,
        "rho_grad_p": rho_grad_p,
        "rho_grad_n_square": rho_grad_n**2,
        "rho_grad_p_square": rho_grad_p**2,
        "rho_grad_square": (rho_grad_n + rho_grad_p)**2,
        "nu_n": 0.1*rho_n,
        "nu_p": 0.1*rho_p,
        "kappa_n": np.full(n, -36630.4),
        "kappa_p": np.full(n, -45207.2),
        "x": np.linspace(-0.5, 10., n),
        # 1D density grid for the functions that differentiate with np.gradient
        "rho_1d": np.linspace(0.01, 0.3, 60),
    }


def _case(target, args, extra=(), scalar=True):
    """One golden case: ``target(*inputs[args], *extra)``.

    Args:
        target (str): ``"module:function"``, module relative to ``libnest``
        args (tuple of str): keys of the input arrays, in call order
        extra (tuple): literal positional arguments appended after the inputs
        scalar (bool): also evaluate with Python-float inputs, point by point
    """
    name = target + ("[" + ",".join(map(str, extra)) + "]" if extra else "")
    return name, {"target": target, "args": list(args), "extra": list(extra),
                  "scalar": scalar}


NP = ("rho_n", "rho_p")
CASES = dict([
    # --- pairing fields
    _case("bsk:neutron_pairing_field", ("rho_n",)),
    _case("bsk:symmetric_pairing_field", NP),
    _case("bsk:neutron_ref_pairing_field", NP),
    _case("bsk:proton_ref_pairing_field", NP),
    _case("bsk:neutron_ref_pairing_field_eq2", NP),
    _case("bsk:proton_ref_pairing_field_eq2", NP),
    _case("bsk:neutron_ref_pairing_field_eq3", NP),
    _case("bsk:proton_ref_pairing_field_eq3", NP),
    _case("bsk:neutron_ref_pairing_field_eq6", NP),
    _case("bsk:proton_ref_pairing_field_eq6", NP),
    # --- thermodynamics
    _case("bsk:energy_per_nucleon", NP),
    _case("bsk:testMe", ("rho",)),
    _case("bsk:energy_per_nucleon_n", ("rho",)),
    _case("bsk:derivative_energy_per_nucleon_n", ("rho",)),
    _case("bsk:second_derivative_energy_per_nucleon_n", ("rho",)),
    _case("bsk:pressure_n", ("rho",)),
    _case("bsk:pressure_derivative_n", ("rho",)),
    _case("bsk:epsilon_derivative_n", ("rho",)),
    _case("bsk:speed_of_sound_n", ("rho",)),
    _case("bsk:numerical_derivative_energy_per_nucleon_n", ("rho_1d",), scalar=False),
    _case("bsk:numerical_second_derivative_energy_per_nucleon_n", ("rho_1d",), scalar=False),
    _case("bsk:numerical_derivative_epsilon_n", ("rho_1d",), scalar=False),
    _case("bsk:numerical_pressure_n", ("rho_1d",), scalar=False),
    _case("bsk:numerical_derivative_pressure_n", ("rho_1d",), scalar=False),
    _case("bsk:numerical_speed_of_sound_n", ("rho_1d",), scalar=False),
    # --- effective masses and mean fields
    _case("bsk:isoscalarM", NP),
    _case("bsk:isovectorM", NP),
    _case("bsk:effMn", NP),
    _case("bsk:effMp", NP),
    _case("bsk:U_q", NP, ("n",)),
    _case("bsk:U_q", NP, ("p",)),
    _case("bsk:B_q", NP, ("n",)),
    _case("bsk:B_q", NP, ("p",)),
    # --- energy functional
    _case("bsk:C0_rho", NP),
    _case("bsk:C1_rho", NP),
    _case("bsk:C0_tau", NP),
    _case("bsk:C1_tau", NP),
    _case("bsk:epsilon_rho_np", NP),
    _case("bsk:epsilon_tau_np", NP + ("tau_n", "tau_p", "jsum2", "jdiff2")),
    _case("bsk:epsilon_delta_rho_np",
          NP + ("rho_grad_n_square", "rho_grad_p_square", "rho_grad_square")),
    _case("bsk:epsilon_pi_np",
          NP + ("rho_grad_n", "rho_grad_p", "nu_n", "nu_p", "kappa_n", "kappa_p")),
    _case("bsk:epsilon_np",
          NP + ("rho_grad_n", "rho_grad_p", "tau_n", "tau_p", "jsum2", "jdiff2",
                "nu_n", "nu_p", "kappa_n", "kappa_p")),
    _case("bsk:v_pi", NP, ("n",)),
    _case("bsk:v_pi", NP, ("p",)),
    _case("bsk:I", NP, ("n",)),
    _case("bsk:I", NP, ("p",)),
    _case("bsk:Lambda", ("x",)),
    # --- functional-dependent helpers in definitions (move to the functional in E2)
    _case("definitions:mu_q", NP, ("n",)),
    _case("definitions:mu_q", NP, ("p",)),
    _case("definitions:xiBCS", ("kF_n",)),
    _case("definitions:E_minigap_rho_n", ("rho_n",)),
])

# Public API that must stay importable with compatible signatures (E2 compat layer).
API_FUNCTIONS = {
    "bsk": None,                    # None = every public function of the module
    "definitions": ["mu_q", "xiBCS", "E_minigap_rho_n"],
}
API_CONSTANTS = {
    "bsk": ["T0", "T1", "T2", "T3", "T4", "T5", "X0", "X1", "T2X2", "X3", "X4", "X5",
            "ALPHA", "BETA", "GAMMA", "YW", "FNP", "FNM", "FPP", "FPM",
            "KAPPAN", "KAPPAP"],
}


def resolve_module(module_name):
    """Returns the module ``libnest.<module_name>``."""
    return importlib.import_module("libnest." + module_name)


def resolve(target):
    """Returns the callable for ``"module:function"`` (module inside ``libnest``)."""
    module, func = target.split(":")
    return getattr(resolve_module(module), func)


def public_functions(module_name):
    """Returns {name: parameter names} of the public functions defined in a module."""
    module = resolve_module(module_name)
    names = API_FUNCTIONS[module_name]
    if names is None:
        names = [n for n, f in inspect.getmembers(module, inspect.isfunction)
                 if f.__module__ == module.__name__ and not n.startswith("_")]
    return {n: list(inspect.signature(getattr(module, n)).parameters) for n in sorted(names)}


def evaluate(case, inputs):
    """Evaluates one case.

    Returns:
        dict: ``array`` (list), ``mask`` (list of bool or None, for masked arrays),
        and, if the case allows scalars, from point-by-point calls with Python
        floats: ``scalar`` (list; None where the call raised), ``scalar_ndim``
        (list of int) and ``scalar_error`` (exception class name or None).
    """
    func = resolve(case["target"])
    args = [np.asarray(inputs[k], dtype=float) for k in case["args"]]
    extra = case["extra"]
    with np.errstate(all="ignore"):
        out = func(*args, *extra)
        result = {"array": np.ma.getdata(out).astype(float).ravel().tolist(),
                  "mask": (np.ma.getmaskarray(out).ravel().tolist()
                           if np.ma.isMaskedArray(out) else None)}
        if case["scalar"]:
            values, ndims, errors = [], [], []
            for i in range(args[0].size):
                try:
                    r = func(*[float(a[i]) for a in args], *extra)
                except Exception as e:  # recorded: the refactor must keep it
                    values.append(None)
                    ndims.append(None)
                    errors.append(type(e).__name__)
                    continue
                values.append(float(r))
                ndims.append(int(np.ndim(r)))
                errors.append(None)
            result["scalar"] = values
            result["scalar_ndim"] = ndims
            result["scalar_error"] = errors
    return result
