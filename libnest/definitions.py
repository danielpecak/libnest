#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Info: who, where, when
# @author Daniel Pęcak <Daniel.Pecak@pw.edu.pl>
# Warsaw Technical University, Université Libre de Bruxelles
# On leave: Institute of Physics, Polish Academy of Sciences, Warsaw
# March 2022, Brussels
"""
Module: Definitions
===================


"""

import numpy as np
import math
# import libnest.bsk
# import libnest.units as units
from libnest import units
from libnest.units import HBARC, DENSEPSILON, NUMZERO
from libnest.units import MN, MP, HBAR2M_n, HBAR2M_p
from libnest.units import DENSEPSILON


def rho2kf(rho):
    """
    Returns wavevector kF based on density rho.

    It uses the relation for a uniform Fermi system and yields:

    .. math::

        k_F = (3 \\pi^2 \\rho )^{1/3}.

    If we have a two-component mixture of protons and neutrons, both
    having two spin components (up and down), then :math:`\\rho` states
    for the density of one isospin, and one spin component only.

    Args:
        rho (float):  density :math:`\\rho` for a single component

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

    If we have a two-component mixture of protons and neutrons, both
    having two spin components (up and down), then :math:`\\rho` states
    for the density of one isospin, and one spin component only.

    Args:
        kF (float):  density :math:`\\rho` for a single component

    Returns:
        float: wavevector :math:`k_F` [fm :sup:`-1`]

    See also:
        :func:`rho2kf`
    """
    return kF**3/(3.*np.pi*np.pi)

def rho2tau(rho):
    """
    Returns kinetic density :math:`\\tau` for uniform Fermi system of density
    :math:`\\rho`.

    .. math::
        \\tau = \\frac{3}{5} \\left(3 \\pi\\right)^{2/3} \\rho^{5/3}

    Args:
        rho (float):  density :math:`\\rho` [fm :sup:`-3`] for a single component

    Returns:
        float: kinetic density :math:`\\tau` [fm :sup:`-3`]
    """
    return 0.6*(3.*np.pi)**(2./3.)*rho**(5./3.)

def rhoEta(rho_n, rho_p):
    """
    Returns total density :math:`\\rho` and difference :math:`\\eta` of
    densities from neutron and proton densities.

    .. math::

        \\rho = \\rho_n + \\rho_p

        \\eta = \\rho_n - \\rho_p

    .. todo::
        Check if I should return rho_n-rho_p OR rho_n-rho_p/(rho_n+rho_p + DENSEPSILON)
        #for now rho_n-rho_p because of the  function: neutron_ref_pairing_field(rho_n, rho_p)

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

    Returns:
        float: pair of total density :math:`\\rho` [fm :sup:`-3`], and density difference
        :math:`\\eta` [fm :sup:`-3`]
    """
    #return rho_n+rho_p, rho_n-rho_p/(rho_n+rho_p + DENSEPSILON)
    return rho_n+rho_p, rho_n-rho_p



def xiBCS(kF, delta=None):
    """
    Calculates the coherence length from BCS theory. If no delta argument is provided,
    it is assumed that we consider pure neutron matter.

    .. math::

        \\xi_\\mathrm{BCS} = \\frac{(\\hbar c)^2 k_F}{\\pi \\Delta Mc^2}

    Args:
        k_F (float): Fermi momentum [fm :sup:`-1`]
        delta (float): pairing field [MeV]

    Returns:
        float: coherence length [fm]
    """
    # ValueError: The truth value of a Series is ambiguous. Use a.empty, a.bool(), a.item(), a.any() or a.all().
    from libnest import bsk
    if delta.all()==None:
        rho   = kf2rho(kF)
        delta = bsk.neutron_pairing_field(rho)
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
        rho_out (float):  density of superlufid neutrons :math:`\\rho_{\\mathrm{out}}` [fm :sup:`-3`]

    Returns:
        float: effective mass :math:`M_{\\mathrm{eff}}` in units neutron mass
    """
    return HBAR2M_n * kF**2

def E_minigap_rho_n(rho_n):
    """
    Returns the energy of minigap :math:`E_{mg}` [MeV] for neutron matter.

    .. todo:: join these two: :func:`E_minigap_delta_n` :func:`E_minigap_rho_n`

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
        :func:`neutron_ref_pairing_field`
        :func:`rho2kf`
        :func:`eF_n`

    """
# https://stackoverflow.com/questions/22095000/how-to-code-a-function-that-accepts-float-list-or-numpy-array
# import numpy as np
#
# def get_lerp_factor(a, x, b):
#     a, x, b = np.asarray(a), np.asarray(x), np.asarray(b)
#     return ((x - a) / (b - a)).clip(0, 1)
    from libnest.bsk import neutron_ref_pairing_field
    delta = neutron_ref_pairing_field(rho_n, 0.)
    return 4./3. * np.abs(delta)**2/eF_n(rho2kf(rho_n))

def E_minigap_delta_n(delta, rho_n):
    """
    Returns the energy of minigap :math:`E_{mg}` [MeV] for neutron matter.

    .. todo:: join these two: :func:`E_minigap_delta_n` :func:`E_minigap_rho_n`


    Parameters
        delta (float): pairing field for neutrons :math:`\\Delta_n` [fm :sup:`-3`]
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

    Returns
        float: energy of minigap :math:`E_{mg}` [MeV]

    See also:
        :func:`eF_n`
    """
    return 4./3. * np.abs(delta)**2/eF_n(rho2kf(rho_n))


def mu_q(rho_n, rho_p, q):
# Eq. taken from S. Goriely, N. Chamel, and J. M. Pearson, Phys. Rev. Lett. 102, 152503 (2009)
    """
    Calculates the chemical potential :math:`\\mu` defined with the wavevector
    :math:`k_F` :cite:`chamel2009pairing`.

    .. math::

        \\mu_q = \\frac{\\hbar^2 k_F^2}{2M_q}

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        q (string): nucleon type choice ('p' - proton, or 'n' - neutron)

    Returns:
         float: chemical potential :math:`\\mu` [MeV]
    """
    if(q=='n'):
        M = effMn(rho_n, rho_p)
        rho = rho_n
    elif(q=='p'):
        M = effMp(rho_n, rho_p)
        rho = rho_p
    else:
        sys.exit('# ERROR: Nucleon q must be either n or p')
    mu_q = HBARC**2*rho2kf(rho)**2/(2.*M)
    i = np.where(mu_q==0)
    mu_q[i] = NUMZERO
    j = np.where(rho==0)
    mu_q[j] = NUMZERO
    return mu_q



if __name__ == '__main__':
    pass
