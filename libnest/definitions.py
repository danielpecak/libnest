# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Daniel Pęcak
"""
Module: Definitions
===================
Handy functions from different fields of physics

List of functions
-----------------
"""

import numpy as np
import math
# import libnest.bsk
# import libnest.units as units
from libnest import units
from libnest.units import DENSEPSILON
from libnest.units import MN, HBAR2M_n


def _functional(functional):
    """Returns the functional to use: BSk31 by default, or one given by name or object."""
    if functional is None or isinstance(functional, str):
        # imported here: libnest.edf imports this module
        from libnest.edf import get_functional
        return get_functional(functional or "BSk31")
    return functional


def rho2kf(rho):
    """
    Returns wavevector kF based on density rho.

    It uses the relation for a uniform Fermi system and yields:

    .. math::

        k_F = (3 \\pi^2 \\rho )^{1/3}.

    Here :math:`\\rho` is the density of a single nucleon species (neutrons
    or protons), summed over both spin components. For symmetric nuclear
    matter of total density :math:`\\rho` pass :math:`\\rho/2`.

    Args:
        rho (float): density :math:`\\rho_q` of one nucleon species, both spin
            components [fm :sup:`-3`]

    Returns:
        float: wavevector :math:`k_F` [fm :sup:`-1`]

    See also:
        :func:`kf2rho`
    """
    return (3.*np.pi*np.pi*rho)**(1./3.)


def kf2rho(kF):
    """
    Returns rho based on wavevector kF.

    It uses the relation for a uniform Fermi system and yields:

    .. math::

        \\rho = \\frac{k_F^3}{3 \\pi^2}.

    The result is the density of a single nucleon species (neutrons or
    protons), summed over both spin components.

    Args:
        kF (float): Fermi wavevector :math:`k_F` [fm :sup:`-1`]

    Returns:
        float: density :math:`\\rho_q` of one nucleon species, both spin
        components [fm :sup:`-3`]

    See also:
        :func:`.rho2kf`
    """
    return kF**3/(3.*np.pi*np.pi)

def rho2tau(rho):
    """
    Returns kinetic density :math:`\\tau` for uniform Fermi system of density
    :math:`\\rho`.

    .. math::
        \\tau = \\frac{3}{5} k_F^2 \\rho = \\frac{3}{5} \\left(3 \\pi^2\\right)^{2/3} \\rho^{5/3}

    Args:
        rho (float): density :math:`\\rho_q` of one nucleon species, both spin
            components [fm :sup:`-3`]

    Returns:
        float: kinetic density :math:`\\tau_q` [fm :sup:`-5`]
    """
    return 0.6*(3.*np.pi*np.pi)**(2./3.)*rho**(5./3.)

def rhoEta(rho_n, rho_p):
    """
    Returns total density :math:`\\rho` and difference :math:`\\eta` of
    densities from neutron and proton densities.

    .. math::

        \\rho = \\rho_n + \\rho_p

        \\eta = \\rho_n - \\rho_p

    .. note::
        The function returns the absolute difference :math:`\\eta = \\rho_n - \\rho_p`,
        not the relative asymmetry :math:`\\eta/(\\rho_n + \\rho_p)`.
        This is consistent with usage in :func:`.neutron_ref_pairing_field` and
        :func:`.proton_ref_pairing_field` where the ratio :math:`\\eta/\\rho` is
        computed explicitly.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

    Returns:
        float: pair of total density :math:`\\rho` [fm :sup:`-3`], and density difference
        :math:`\\eta` [fm :sup:`-3`]
    """
    #return rho_n+rho_p, rho_n-rho_p/(rho_n+rho_p + DENSEPSILON)
    return rho_n+rho_p, rho_n-rho_p



def xiBCS(kF, delta=None, functional=None):
    """
    Calculates the coherence length from BCS theory. If no delta argument is provided,
    it is assumed that we consider pure neutron matter, with the neutron pairing field
    of the functional.

    .. math::

        \\xi_\\mathrm{BCS} = \\frac{(\\hbar c)^2 k_F}{\\pi \\Delta Mc^2}

    Args:
        k_F (float): Fermi momentum [fm :sup:`-1`]
        delta (float): pairing field [MeV]
        functional (SkyrmeFunctional or str): functional for the default pairing
            field; BSk31 if not given

    Returns:
        float: coherence length [fm]
    """
    # ValueError: The truth value of a Series is ambiguous. Use a.empty, a.bool(), a.item(), a.any() or a.all().
    kF = np.asarray(kF, dtype=float)
    if delta is None:
        rho   = kf2rho(kF)
        delta = _functional(functional).neutron_pairing_field(rho)
    return units.HBARC**2*kF/(np.pi*delta*units.MN)


def vsf(r):
    """
    Calculates the velocity based on gradient of the pairing field gradient.

    .. math::

        v_\\mathrm{sf} = \\frac{\\hbar c}{M} \\frac{1}{2r}

    Args:
        r (float): distance from the center of a vortex in femtometers [fm]

    Returns:
        float: velocity in units of percentage of speed of light [c]
    """
    return units.ALPHA/units.MN*0.5/abs(r)

def vLandau(delta,kF):
    """
    Returns Landau velocity.

    Landau velocity shows at which velocity the superfluid medium starts to be
    excited: phonons appear.

    .. math::

        v_L = \\frac{\\Delta}{\\hbar k_F} c

    Args:
        delta (float): :math:`\\Delta` [MeV]
        kF (float):  :math:`k_F` [MeV]

    Returns:
        float: Landau velocity :math:`v_L` in units of speed of light [c]
    """
    return delta/kF/units.HBARC

def vcritical(delta,kF):
    """
    Returns critical velocity for superfluid.

    At this velocity the system is no longer superfluid: the Cooper pairs break.

    .. math::
        v_L = \\frac{e}{2} \\frac{\\Delta}{\\hbar k_F} c

    Args:
        delta (float): :math:`\\Delta` [MeV]
        kF (float):  :math:`k_F` [MeV]

    Returns:
        float: Landau velocity :math:`v_c` in units of speed of light [c]
    """
    return math.e/2*delta/kF/units.HBARC

def superfluidFraction(j,rho,vsf):
    """
    Returns superfluid fraction: how much of the matter is superfluid.

    .. math::

        \\eta_\\mathrm{sf} = \\frac{\\hbar c}{M} \\frac{{ j}}{\\rho} \\frac{1}{v_\\mathrm{sf}}

    Args:
        rho (float): density :math:`\\rho` [fm :sup:`-3`]
        j (float): a three-component current vector :math:`\\vec j` [fm :sup:`-4`]
        vsf (float):  :math:`v_{\\mathrm{sf}}` [c]

    Returns:
        float: superfluid fraction: a number between 0 and 1 [1]
    """
    return units.ALPHA/units.MN*j/(rho*vsf)

def vsf_NV(B,vsf,A):
    """
    Returns the velocity based on the gradient of the pairing field phase.
    However it is adjusted to the entrainment effects (definition by Nicolas Chamel
    Valentin Allard).

    .. math::

        v_\\mathrm{sf}^{NV} = \\frac{\\hbar^2}{2 M B}v_\\mathrm{sf} + \\frac{{ A}}{M}

    Args:
        B (float): mean field potential coming from kinetic energy variation B [MeV fm :sup:`2`]
        vsf (float):  :math:`v_{\\mathrm{sf}}` [c]
        A (float): mean field potential coming from current variation A [MeV fm]

    Returns:
        float: velocity :math:`v_\\mathrm{SF}^{NV}` in units of percentage of speed of light [c]
    """
    return units.hbar22M0/B*vsf + units.ALPHA/units.MN*A/units.HBARC

def v_NV(B,j,rho,A):
    """
    Returns the velocity (mass velocity). It is adjusted to the entrainment
    effects (definition by Nicolas Chamel Valentin Allard).

    .. math::

        v^{\\mathrm{NV}} = \\hbar c \\frac{\\hbar^2}{2 M B} \\frac{{ j}}{\\rho}  + \\frac{{ A}}{M}

    Args:
        B (float): mean field potential coming from kinetic energy variation B [MeV fm :sup:`2`]
        j (float): a three-component current vector :math:`\\vec j` [fm :sup:`-4`]
        rho (float): density :math:`\\rho` [fm :sup:`-3`]
        A (float): mean field potential coming from current variation A [MeV fm]

    Returns:
        float: mass velocity :math:`v_{\\mathrm{NV}}` in units of percentage of speed of light [c]
    """
    return units.ALPHA/units.MN*(units.hbar22M0/B*j/rho + A/units.HBARC)

def eF_n(kF):
    """
    Returns Fermi energy for neutrons based on wavevector kF:

    .. math::

        \\epsilon_F = \\frac{\\hbar^2 k_F^2}{2 M_n},

    where :math:`M_n` is mass of a neutron.

    Args:
        kF (float):  wavevector :math:`k_F`

    Returns:
        float: Fermi energy :math:`\\epsilon_F` [MeV]
    """
    return HBAR2M_n * kF**2

def Meff_hydro(rho_in, rho_out, R):
    """
    Returns effective mass of a spherical nucleus of radius :math:`R` of density :math:`\\rho_{\\mathrm{in}}` immersed in superfluid neutrons of density :math:`\\rho_{\\mathrm{out}}`.
    The formula is from: :cite:`magierski2004medium`

    .. math::

        M_{\\mathrm{eff}} = \\frac{4}{3} \\pi R^3 m_n \\frac{(\\rho_{\\mathrm{in}}-\\rho_{\\mathrm{out}})^2}{\\rho_{\\mathrm{in}}+2\\rho_{\\mathrm{out}}},

    where :math:`m_n` is mass of a neutron.

    Args:
        rho_in (float):  density of nucleus :math:`\\rho_{\\mathrm{in}}` [fm :sup:`-3`]
        rho_out (float):  density of superfluid neutrons :math:`\\rho_{\\mathrm{out}}` [fm :sup:`-3`]

    Returns:
        float: effective mass :math:`M_{\\mathrm{eff}}` in units neutron mass
    """
    return (4./3.)*np.pi*R**3*MN*(rho_in-rho_out)**2/(rho_in+2.*rho_out+DENSEPSILON)

def E_minigap_delta_n(delta, rho_n):
    """
    Returns the energy of minigap :math:`E_{mg}` [MeV] for neutron matter
    given the pairing field and density.

    The minigap energy can be approximated:

    .. math::

        E_\\mathrm{mg} = \\frac{4}{3} \\frac{|\\Delta|^2}{\\varepsilon_F},

    where :math:`\\Delta` is the pairing gap in the system, and :math:`\\varepsilon_F`
    is the Fermi energy.

    Args:
        delta (float): pairing field for neutrons :math:`\\Delta_n` [MeV]
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both
            spin components

    Returns:
        float: energy of minigap :math:`E_{mg}` [MeV]

    See also:
        :func:`.E_minigap_rho_n`
        :func:`.eF_n`
    """
    return 4./3. * np.abs(delta)**2/eF_n(rho2kf(rho_n))

def E_minigap_rho_n(rho_n, functional=None):
    """
    Returns the energy of minigap :math:`E_{mg}` [MeV] for neutron matter,
    calculating the pairing field automatically from density.

    This is a convenience wrapper around :func:`.E_minigap_delta_n` that
    computes the neutron pairing field from the functional (BSk31 by default).

    The minigap energy can be approximated:

    .. math::

        E_\\mathrm{mg} = \\frac{4}{3} \\frac{|\\Delta|^2}{\\varepsilon_F},

    where :math:`\\Delta` is the pairing gap in the system, and :math:`\\varepsilon_F`
    is the Fermi energy.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both
            spin components

    Returns:
        float: energy of minigap :math:`E_{mg}` [MeV]

    See also:
        :func:`.E_minigap_delta_n`
        :func:`.neutron_ref_pairing_field`
        :func:`.rho2kf`
        :func:`.eF_n`

    """
    delta = _functional(functional).neutron_ref_pairing_field(rho_n, 0.)
    return E_minigap_delta_n(delta, rho_n)



def mu_q(rho_n, rho_p, q, functional=None):
    # Eq. taken from S. Goriely, N. Chamel, and J. M. Pearson, Phys. Rev. Lett. 102, 152503 (2009)
    """
    Calculates the chemical potential :math:`\\mu_q` defined with the wavevector
    :math:`k_{F,q}` and the effective mass :math:`M^*_q` :cite:`chamel2009pairing`.

    .. math::

        \\mu_q = \\frac{\\hbar^2 k_{F,q}^2}{2M^*_q} = B_q k_{F,q}^2

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        q (str): nucleon type choice ('p' - proton, or 'n' - neutron)
        functional (SkyrmeFunctional or str): functional that provides the
            effective mass; BSk31 if not given

    Returns:
         float: chemical potential :math:`\\mu` [MeV]
    """
    return _functional(functional).mu_q(rho_n, rho_p, q)


if __name__ == '__main__':
    pass
