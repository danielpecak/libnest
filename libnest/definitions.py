#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Info: who, where, when
# @author Daniel Pęcak <Daniel.Pecak@pw.edu.pl>
# Warsaw Technical University, Université Libre de Bruxelles
# On leave: Institute of Physics, Polish Academy of Sciences, Warsaw
# March 2022, Brussels
# =========== Description
# Script for converting the density units for nuclear matter
# in the ranges typical for neutron stars.
# =========== Usage example
# $ ./units.py
import numpy as np
import math
import libnest.units

def kf(rho):
    """Returns the Fermi momentum kF [MeV] for a density rho [fm^{-3}]"""
    return np.power(3.*math.pi*math.pi*rho,1./3.)

def vsf(r):
    """Returns the velocity based on gradient of the pairing field gradient.
In units of percentage of speed of light [c]."""
    return libnest.units.VUNIT*0.5/abs(r)

def vLandau(delta,kF):
    """Returns Landau velocity."""
    return delta/kF

def vcritical(delta,kF):
    """Returns critical velocity for superfluid."""
    return delta/kF*math.e/2.

def superfluidFraction(j,rho,vsf):
    """Returns superfluid fraction: how much of the matter is superfluid."""
    return libnest.units.VUNIT*j/(rho*vsf)

def vsf_NV(B,vsf,A):
    """Returns the velocity based on the gradient of the pairing field phase.
However it is adjusted to the entrainment effects (definition by Nicolas Chamel
Valentin Allard)."""
    return libnest.units.M0/B*vsf + libnest.units.VUNIT*A/libnest.units.HBARC

def v_NV(B,j,rho,A):
    """Returns the velocity (mass velocity). It is adjusted to the entrainment
effects (definition by Nicolas Chamel Valentin Allard)."""
    return libnest.units.VUNIT*(libnest.units.M0/B*j/rho + A/libnest.units.HBARC)

if __name__ == '__main__':
    pass
