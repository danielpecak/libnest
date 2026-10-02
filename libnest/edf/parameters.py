# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Daniel Pęcak
"""
Parameter sets of Skyrme-type energy density functionals.

A parametrization is data: :class:`SkyrmeParameters` holds the Skyrme part and
:class:`PairingParameters` the pairing part. Both are immutable; derive a modified
set with :func:`dataclasses.replace`, e.g. ``replace(params, t0=-2300.)``.
"""
from dataclasses import dataclass
from typing import Optional

from libnest.units import HBAR2M_n, HBAR2M_p


@dataclass(frozen=True)
class SkyrmeParameters:
    """
    Parameters of an (extended) Skyrme functional.

    Standard Skyrme forces have no density-dependent gradient terms: leave
    ``t4 = t5 = 0`` and ``beta = gamma = 0`` at their defaults.

    Attributes:
        name (str): name of the parametrization, e.g. ``"BSk31"``
        t0 (float): :math:`t_0` [MeV fm :sup:`3`]
        t1 (float): :math:`t_1` [MeV fm :sup:`5`]
        t2 (float): :math:`t_2` [MeV fm :sup:`5`]
        t3 (float): :math:`t_3` [MeV fm :sup:`3+3\\alpha`]
        x0 (float): :math:`x_0` [1]
        x1 (float): :math:`x_1` [1]
        t2x2 (float): the product :math:`t_2 x_2` [MeV fm :sup:`5`]; stored directly
            because BSk forces have :math:`t_2 = 0` with finite :math:`t_2 x_2`
        x3 (float): :math:`x_3` [1]
        alpha (float): density exponent :math:`\\alpha` of the :math:`t_3` term [1]
        t4 (float): :math:`t_4` [MeV fm :sup:`5+3\\beta`]
        t5 (float): :math:`t_5` [MeV fm :sup:`5+3\\gamma`]
        x4 (float): :math:`x_4` [1]
        x5 (float): :math:`x_5` [1]
        beta (float): density exponent :math:`\\beta` of the :math:`t_4` term [1]
        gamma (float): density exponent :math:`\\gamma` of the :math:`t_5` term [1]
        hbar2m_n (float): :math:`\\hbar^2/(2 m_n)` [MeV fm :sup:`2`]
        hbar2m_p (float): :math:`\\hbar^2/(2 m_p)` [MeV fm :sup:`2`]
        yw (float or None): ``YW`` of the original ``bsk`` module; not used in
            uniform-matter formulas
        reference (str): publication that defines the parameters
        doi (str): DOI of that publication
        notes (str): free-form remarks
    """
    name: str
    t0: float
    t1: float
    t2: float
    t3: float
    x0: float
    x1: float
    t2x2: float
    x3: float
    alpha: float
    t4: float = 0.0
    t5: float = 0.0
    x4: float = 0.0
    x5: float = 0.0
    beta: float = 0.0
    gamma: float = 0.0
    hbar2m_n: float = HBAR2M_n
    hbar2m_p: float = HBAR2M_p
    yw: Optional[float] = None
    reference: str = ""
    doi: str = ""
    notes: str = ""


@dataclass(frozen=True)
class PairingParameters:
    """
    Parameters of the pairing part of a functional.

    Attributes:
        cutoff (float): single-particle energy cutoff :math:`\\varepsilon_\\Lambda`
            [MeV] used in the pairing strength :math:`v^{\\pi}_q`
        kappa_n (float): gradient-pairing coefficient :math:`\\kappa_n` [MeV fm :sup:`8`]
        kappa_p (float): gradient-pairing coefficient :math:`\\kappa_p` [MeV fm :sup:`8`]
        f_np (float): pairing renormalization factor :math:`f^+_n` [1]
        f_nm (float): pairing renormalization factor :math:`f^-_n` [1]
        f_pp (float): pairing renormalization factor :math:`f^+_p` [1]
        f_pm (float): pairing renormalization factor :math:`f^-_p` [1]
    """
    cutoff: float = 6.5
    kappa_n: float = 0.0
    kappa_p: float = 0.0
    f_np: float = 1.0
    f_nm: float = 1.0
    f_pp: float = 1.0
    f_pm: float = 1.0
