# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Daniel Pęcak
"""
Reference pairing gaps of uniform matter.

The pairing strength :math:`v^{\\pi q}` of a functional (see
:meth:`.SkyrmeFunctional.v_pi`) is fixed by requiring that it reproduces, in uniform
matter, a set of reference :math:`{}^1S_0` gaps :math:`\\Delta_q(\\rho_n,\\rho_p)`.
They are built from two ingredients, both selected by name in
:class:`~libnest.edf.PairingParameters`:

* a **gap model** (``gap_model``, see :data:`GAP_MODELS`): the gaps in pure neutron
  matter, :math:`\\Delta_{\\mathrm{NeuM}}(k_F)`, and in symmetric matter,
  :math:`\\Delta_{\\mathrm{SM}}(k_F)`;
* an **interpolation scheme** (``scheme``, see :data:`SCHEMES`) that combines them into
  :math:`\\Delta_n` and :math:`\\Delta_p` at arbitrary asymmetry.

All schemes are isospin symmetric, :math:`\\Delta_p(\\rho_n,\\rho_p) =
\\Delta_n(\\rho_p,\\rho_n)`, so a scheme is written once, for species :math:`q`, with
:math:`\\bar q` the other species and the asymmetry

.. math::

    \\delta_q = \\frac{\\rho_q - \\rho_{\\bar q}}{\\rho}, \\qquad \\rho = \\rho_n + \\rho_p .

Every scheme reduces to :math:`\\Delta_{\\mathrm{SM}}` in symmetric matter and to
:math:`\\Delta_{\\mathrm{NeuM}}(\\rho_n)` for neutrons in pure neutron matter.
"""
from dataclasses import dataclass
from typing import Callable

import numpy as np

from libnest.units import DENSEPSILON, NUMZERO
from libnest.definitions import rho2kf, rhoEta


def _clamp_nonnegative(pairing):
    """Floors an interpolated pairing field at (numerical) zero.

    The SM/NeuM interpolations carry a minus sign in one of the two terms, so
    for sufficiently asymmetric matter the result can dip below zero.

    Args:
        pairing (float or numpy.ndarray): raw interpolated pairing field [MeV]

    Returns:
        float or numpy.ndarray: the same, with negative values replaced by NUMZERO
    """
    pairing = np.maximum(np.asarray(pairing, dtype=float), NUMZERO)
    if pairing.shape == ():
        return float(pairing)
    return pairing


def _close_above(delta, kF, kF_max):
    """Sets the gap to (numerical) zero where ``kF > kF_max``.

    Scalar input returns early: ``np.where`` on a 0-d array raises in NumPy 2.
    """
    if delta.shape == ():
        if kF > kF_max:
            return float(NUMZERO)
        return float(delta)
    delta[np.where(kF > kF_max)] = NUMZERO
    return delta


# ================================
#   Gap models
# ================================
def cao2006_neutron_matter(kF):
    r"""
    Gap in pure neutron matter: fit (NeST.pdf Eq. 5.12; arXiv:2411.08007, Table 5) to
    the EBHF gaps of Cao, Lombardo, Schuck, PRC 74, 064301 (2006), with polarization and
    self-energy. For :math:`k_F > 1.38` fm :sup:`-1` it returns (numerical) zero.

    .. math::

        \Delta_{\mathrm{NeuM}}(k_F) = \frac{3.37968 k_F^2}{k_F^2+0.556092^2}
        \frac{(k_F-1.38236)^2}{(k_F-1.38236)^2+0.327517^2}

    Args:
        kF (float or numpy.ndarray): neutron Fermi wavevector :math:`k_F` [fm :sup:`-1`]

    Returns:
        float or numpy.ndarray: :math:`\Delta_{\mathrm{NeuM}}` [MeV]
    """
    kF = np.asarray(kF, dtype=float)
    delta = 3.37968*(kF**2)*((kF-1.38236)**2)/(((kF**2)+(0.556092**2))*
                                              ((kF-1.38236)**2+(0.327517**2)))
    return _close_above(delta, kF, 1.38)


def cao2006_symmetric_matter(kF):
    r"""
    Gap in symmetric matter: fit (NeST.pdf Eq. 5.11; arXiv:2411.08007, Table 5) to the
    EBHF gaps of Cao, Lombardo, Schuck, PRC 74, 064301 (2006), with polarization and
    self-energy. For :math:`k_F > 1.31` fm :sup:`-1` it returns (numerical) zero.

    .. math::

        \Delta_{\mathrm{SM}}(k_F) = \frac{11.5586 k_F^2}{k_F^2 + 0.489932^2}
        \frac{(k_F - 1.3142)^2}{(k_F - 1.3142)^2 + 0.906146^2}

    Args:
        kF (float or numpy.ndarray): Fermi wavevector of one species,
            :math:`k_F = (3\pi^2\rho/2)^{1/3}` [fm :sup:`-1`]

    Returns:
        float or numpy.ndarray: :math:`\Delta_{\mathrm{SM}}` [MeV]
    """
    kF = np.asarray(kF, dtype=float)
    delta = 11.5586*(kF**2)*((kF-1.3142)**2)/(((kF**2)+(0.489932**2))*
                                             (((kF-1.3142)**2)+(0.906146**2)))
    return _close_above(delta, kF, 1.31)


@dataclass(frozen=True)
class GapModel:
    """
    Reference gaps in pure neutron matter and in symmetric matter.

    Attributes:
        neutron_matter_kf (callable): :math:`\\Delta_{\\mathrm{NeuM}}(k_F)` [MeV]
        symmetric_matter_kf (callable): :math:`\\Delta_{\\mathrm{SM}}(k_F)` [MeV], with
            :math:`k_F` of one species
        reference (str): where the gaps come from
    """
    neutron_matter_kf: Callable
    symmetric_matter_kf: Callable
    reference: str = ""

    def neutron_matter(self, rho_q):
        """
        Args:
            rho_q (float or numpy.ndarray): density :math:`\\rho_q` [fm :sup:`-3`] of
                one species, both spin components

        Returns:
            float or numpy.ndarray: :math:`\\Delta_{\\mathrm{NeuM}}` at
            :math:`k_F = (3\\pi^2\\rho_q)^{1/3}` [MeV]
        """
        return self.neutron_matter_kf(rho2kf(rho_q))

    def symmetric_matter(self, rho_n, rho_p):
        """
        Args:
            rho_n (float or numpy.ndarray): neutron density :math:`\\rho_n` [fm :sup:`-3`]
            rho_p (float or numpy.ndarray): proton density :math:`\\rho_p` [fm :sup:`-3`]

        Returns:
            float or numpy.ndarray: :math:`\\Delta_{\\mathrm{SM}}` at the per-species
            :math:`k_F = (3\\pi^2\\rho/2)^{1/3}`, :math:`\\rho = \\rho_n + \\rho_p` [MeV]
        """
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        # per-species Fermi momentum: in symmetric matter rho_n = rho_p = rho/2
        # (Chamel, Goriely, Pearson, PRC 80, 065804 (2009))
        return self.symmetric_matter_kf(rho2kf(0.5*(rho_n+rho_p)))


#: Gap models by name (:class:`GapModel`).
GAP_MODELS = {
    "cao2006": GapModel(
        cao2006_neutron_matter, cao2006_symmetric_matter,
        reference="L. G. Cao, U. Lombardo, P. Schuck, Phys. Rev. C 74, 064301 (2006); "
                  "fits for BSk30-32, listed in arXiv:2411.08007, Table 5"),
}


# ================================
#   Interpolation schemes
# ================================
def chamel2009(gaps, rho_q, rho_qbar):
    r"""
    Linear interpolation of Chamel, Goriely, Pearson, PRL 102, 152503 (2009) (NeST.pdf
    Eq. 5.10; arXiv:2411.08007, Eq. 2), used by the BSk family since BSk17:

    .. math::

        \Delta_q = \Delta_{\mathrm{SM}}(\rho)\,(1 - |\delta_q|)
        + \delta_q \frac{\rho_q}{\rho} \Delta_{\mathrm{NeuM}}(\rho_q) .

    For :math:`\delta_q` close to :math:`-1` the formula turns negative; the result is
    floored at (numerical) zero.

    Args:
        gaps (GapModel): reference gaps of SM and NeuM
        rho_q (float or numpy.ndarray): density of species :math:`q` [fm :sup:`-3`]
        rho_qbar (float or numpy.ndarray): density of the other species [fm :sup:`-3`]

    Returns:
        float or numpy.ndarray: :math:`\Delta_q` [MeV]
    """
    rho_q = np.asarray(rho_q, dtype=float)
    rho_qbar = np.asarray(rho_qbar, dtype=float)
    rho, eta = rhoEta(rho_q, rho_qbar)
    rho = np.asarray(rho + DENSEPSILON, dtype=float)
    return _clamp_nonnegative(gaps.symmetric_matter(rho_q, rho_qbar)*(1-np.abs(eta/rho))
                              + gaps.neutron_matter(rho_q)*rho_q/rho*eta/rho)


def bskg3_eq3(gaps, rho_q, rho_qbar):
    r"""
    Linear interpolation in :math:`|\delta|` of BSkG3 (Grams et al., EPJA 59, 270
    (2023), Eq. 14; arXiv:2411.08007, Eq. 3); both terms are non-negative:

    .. math::

        \Delta_q = \Delta_{\mathrm{SM}}(\rho)\,(1 - |\delta_q|)
        + |\delta_q|\, \Delta_{\mathrm{NeuM}}(\rho_q) .

    Args:
        gaps (GapModel): reference gaps of SM and NeuM
        rho_q (float or numpy.ndarray): density of species :math:`q` [fm :sup:`-3`]
        rho_qbar (float or numpy.ndarray): density of the other species [fm :sup:`-3`]

    Returns:
        float or numpy.ndarray: :math:`\Delta_q` [MeV]
    """
    rho_q = np.asarray(rho_q, dtype=float)
    rho_qbar = np.asarray(rho_qbar, dtype=float)
    rho, eta = rhoEta(rho_q, rho_qbar)
    rho = np.asarray(rho + DENSEPSILON, dtype=float)
    delta = eta/rho
    return _clamp_nonnegative((1-np.abs(delta))*gaps.symmetric_matter(rho_q, rho_qbar)
                              + np.abs(delta)*gaps.neutron_matter(rho_q))


def bskg4_eq6(gaps, rho_q, rho_qbar):
    r"""
    Geometric interpolation of BSkG4 (Grams et al., arXiv:2411.08007, Eq. 6), which
    follows the isospin splitting of the effective mass:

    .. math::

        \Delta_q = \Delta_{\mathrm{NeuM}}(\rho_q)
        \left[ \frac{\Delta_{\mathrm{SM}}(\rho)}{\Delta_{\mathrm{NeuM}}(\rho/2)}
        \right]^{1 - \delta_q} .

    Where both reference gaps are closed the ratio is that of the numerical zeros
    (1); the C code floors both at :math:`10^{-8}` MeV instead, which changes the result
    by less than :math:`10^{-6}` MeV, only where the gap is closed (``TODO_edf.md`` §5).

    Args:
        gaps (GapModel): reference gaps of SM and NeuM
        rho_q (float or numpy.ndarray): density of species :math:`q` [fm :sup:`-3`]
        rho_qbar (float or numpy.ndarray): density of the other species [fm :sup:`-3`]

    Returns:
        float or numpy.ndarray: :math:`\Delta_q` [MeV]
    """
    rho_q = np.asarray(rho_q, dtype=float)
    rho_qbar = np.asarray(rho_qbar, dtype=float)
    rho, eta = rhoEta(rho_q, rho_qbar)
    rho = np.asarray(rho + DENSEPSILON, dtype=float)
    delta = eta/rho
    return _clamp_nonnegative(gaps.neutron_matter(rho_q)*(
        (gaps.symmetric_matter(rho_q, rho_qbar)/gaps.neutron_matter(rho/2))**(1-delta)))


#: Interpolation schemes by name; each is called as ``scheme(gaps, rho_q, rho_qbar)``.
SCHEMES = {
    "chamel2009": chamel2009,
    "bskg3_eq3": bskg3_eq3,
    "bskg4_eq6": bskg4_eq6,
}
