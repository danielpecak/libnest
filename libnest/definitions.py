#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Info: who, where, when
# @author Daniel Pęcak <Daniel.Pecak@pw.edu.pl>
# Warsaw Technical University, Université Libre de Bruxelles
# On leave: Institute of Physics, Polish Academy of Sciences, Warsaw
# March 2022, Brussels
"""
definitions.py
==============
"""

import numpy as np
import math
import libnest.units as units
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

    .. todo::
        formula

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

    .. todo::
        Check if I should return rho_n-rho_p OR rho_n-rho_p/(rho_n+rho_p + DENSEPSILON)
        > na razie oba chyba działają tak samo

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

    Returns:
        float: pair of total density :math:`\\rho` [fm :sup:`-3`], and density difference
        :math:`\\eta` [fm :sup:`-3`]
    """
    return rho_n+rho_p, rho_n-rho_p/(rho_n+rho_p + DENSEPSILON)

def vsf(r):
    """
    Calculates the velocity based on gradient of the pairing field gradient.

    Args:
        r (float): distance from the center of a vortex in femtometers [fm]

    Returns:
        float: velocity in units of percentage of speed of light [c]
    """
    return units.VUNIT*0.5/abs(r)

def vLandau(delta,kF):
    """
    Returns Landau velocity.

    Landau velocity shows at which velocity the superfluid medium starts to be
    excited.

    .. math::
        v_L = \\frac{\\Delta}{\\hbar k_F} c

    Args:
        delta (float): :math:`\\Delta` [MeV]
        kF (float):  :math:`k_F` [MeV]

    Returns:
        float: Landay velocity :math:`v_L` in units of speed of light [c]
    """
    return delta/kF

def vcritical(delta,kF):
    """
    Returns critical velocity for superfluid.

    At this velocity the system is no longer superfluid.

    .. math::
        v_L = e \\frac{\\Delta}{\\hbar k_F} c

    Args:
        delta (float): :math:`\\Delta` [MeV]
        kF (float):  :math:`k_F` [MeV]

    Returns:
        float: Landau velocity :math:`v_c` in units of speed of light [c]
    """
    return delta/kF*math.e/2.

def superfluidFraction(j,rho,vsf):
    """
    Returns superfluid fraction: how much of the matter is superfluid.

    Args:
        rho (float): density :math:`\\rho` [fm :sup:`-3`]
        j (float): a three-component current vector :math:`\\vec j` [fm :sup:`-4`]
        vsf (float):  :math:`v_{\\mathrm{sf}}` [c]

    Returns:
        float: superfluid fraction: a number between 0 and 1 [1]
    """
    return units.VUNIT*j/(rho*vsf)

def vsf_NV(B,vsf,A):
    """
    Returns the velocity based on the gradient of the pairing field phase.
    However it is adjusted to the entrainment effects (definition by Nicolas Chamel
    Valentin Allard).

    Args:
        B (float): mean field potential coming from kinetic energy variation B [MeV fm :sup:`2`]
        vsf (float):  :math:`v_{\\mathrm{sf}}` [c]
        A (float): mean field potential coming from current variation A [MeV fm]

    Returns:
        float: velocity :math:`v_{\\mathrm{SF NV}}` in units of percentage of speed of light [c]
    """
    return units.hbar22M0/B*vsf + units.VUNIT*A/units.HBARC

def v_NV(B,j,rho,A):
    """
    Returns the velocity (mass velocity). It is adjusted to the entrainment
    effects (definition by Nicolas Chamel Valentin Allard).

    Args:
        B (float): mean field potential coming from kinetic energy variation B [MeV fm :sup:`2`]
        j (float): a three-component current vector :math:`\\vec j` [fm :sup:`-4`]
        rho (float): density :math:`\\rho` [fm :sup:`-3`]
        A (float): mean field potential coming from current variation A [MeV fm]

    Returns:
        float: mass velocity :math:`v_{\\mathrm{NV}}` in units of percentage of speed of light [c]     
    """
    return units.VUNIT*(units.hbar22M0/B*j/rho + A/units.HBARC)


if __name__ == '__main__':
    pass
