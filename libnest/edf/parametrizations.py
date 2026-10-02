# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Daniel Pęcak
"""
Registry of the available functionals.

Parameter sets live in their family modules (:mod:`libnest.bsk` for the BSk family,
:mod:`libnest.bskg` for the BSkG family), which call :func:`register` when they are
imported. :func:`get_functional` imports the family modules on first use.
"""
import importlib

from libnest.edf.skyrme import SkyrmeFunctional

# Family modules that register their parameter sets on import.
FAMILY_MODULES = ("libnest.bsk", "libnest.bskg")

_REGISTRY = {}


def register(params, pairing=None):
    """
    Creates a functional from parameter sets and makes it available by name.

    Args:
        params (SkyrmeParameters): Skyrme part; ``params.name`` is the lookup key
            (case-insensitive)
        pairing (PairingParameters): pairing part; defaults to ``PairingParameters()``

    Returns:
        SkyrmeFunctional: the registered functional

    Raises:
        ValueError: if a different functional is already registered under that name
    """
    functional = SkyrmeFunctional(params, pairing)
    key = params.name.lower()
    old = _REGISTRY.get(key)
    if old is not None and (old.params, old.pairing) != (functional.params, functional.pairing):
        raise ValueError("a different functional named {!r} is already registered"
                         .format(params.name))
    _REGISTRY.setdefault(key, functional)
    return _REGISTRY[key]


def _load_families():
    for module in FAMILY_MODULES:
        importlib.import_module(module)


def get_functional(name):
    """
    Returns a registered functional by name.

    Args:
        name (str): e.g. ``"BSk31"``; case-insensitive

    Returns:
        SkyrmeFunctional: the functional

    Raises:
        KeyError: if no functional of that name is registered
    """
    _load_families()
    try:
        return _REGISTRY[name.lower()]
    except KeyError:
        raise KeyError("unknown functional {!r}; available: {}".format(
            name, ", ".join(available_functionals()))) from None


def available_functionals():
    """
    Returns the names of all registered functionals.

    Returns:
        list of str: names, sorted case-insensitively
    """
    _load_families()
    return sorted((f.name for f in _REGISTRY.values()), key=str.lower)
