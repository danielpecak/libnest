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


def rho2kf(rho):
    """Returns wavevector kF based on density rho.

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
    """Returns rho based on wavevector kF.

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
    """Returns kinetic density tau for uniform Fermi system of density rho."""
    return 0.6*(3.*np.pi)**(2./3.)*rho**(5./3.)

def rhoEta(rho_n, rho_p):
    """Returns total density and difference of densities from neutron and
    proton densities."""
    return rho_n+rho_p, rho_n-rho_p


def vsf(r):
    """Calculates the velocity based on gradient of the pairing field gradient.

    Args:
        r (float): distance from the center of a vortex in femtometers [fm]

    Returns:
        float: velocity in units of percentage of speed of light [c]
    """
    return units.VUNIT*0.5/abs(r)

def vLandau(delta,kF):
    """Returns Landau velocity."""
    return delta/kF

def vcritical(delta,kF):
    """Returns critical velocity for superfluid."""
    return delta/kF*math.e/2.

def superfluidFraction(j,rho,vsf):
    """Returns superfluid fraction: how much of the matter is superfluid."""
    return units.VUNIT*j/(rho*vsf)

def vsf_NV(B,vsf,A):
    """Returns the velocity based on the gradient of the pairing field phase.
However it is adjusted to the entrainment effects (definition by Nicolas Chamel
Valentin Allard)."""
    return units.hbar22M0/B*vsf + units.VUNIT*A/units.HBARC

def v_NV(B,j,rho,A):
    """Returns the velocity (mass velocity). It is adjusted to the entrainment
effects (definition by Nicolas Chamel Valentin Allard)."""
    return units.VUNIT*(units.hbar22M0/B*j/rho + A/units.HBARC)




#def isovector_effective_mass(M_q, rho, rho_q):
#    return


if __name__ == '__main__':
    pass
