# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Daniel Pęcak
"""
Energy density functionals: parameter sets, one implementation of the formulas,
and a registry to pick a functional by name.

Example:
    >>> from libnest.edf import get_functional
    >>> f = get_functional("BSk31")
    >>> e = f.energy_per_nucleon(0.08, 0.08)   # symmetric matter, E/A [MeV]
"""
from libnest.edf.parameters import SkyrmeParameters, PairingParameters
from libnest.edf.skyrme import SkyrmeFunctional
from libnest.edf.parametrizations import register, get_functional, available_functionals

__all__ = ["SkyrmeParameters", "PairingParameters", "SkyrmeFunctional",
           "register", "get_functional", "available_functionals"]
