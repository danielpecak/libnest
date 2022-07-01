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

def F_x(x, eta):
    """Returns the function Fx(eta) later used in the calculation of energy per
    nucleon of infinite nuclear matter of density rho and asymmetry eta"""
    return 0.5*((1+eta)**x+(1-eta)**x)

def assymetry_eta(rho, rho_n, rho_p):
    """Returns asymmetry between neutrons and protons eta for matter of
    density rho"""
    return (rho_n-rho_p)/rho

def energy_per_nucleon(kF, rho, F_x, eta):
    """Returns the energy per nucleon on infinite nuclear matter of density rho
    and asymmetry eta, in MeV."""
    return (3*libnest.units.HBARC**2/20*kF**2*((np.power(1+eta,5/3)/libnest.units.MN
                                                +np.power(1-eta,5/3)
                                                          /libnest.units.MP))
            + libnest.units.T0/8*rho*(3-(2*libnest.units.X0 + 1)*eta**2) 
            + 3*libnest.units.T1/40*rho*kF**2*((2+libnest.units.X1)*F_x(5/3,eta)
                                               -(1/2+libnest.units.X1)*F_x(8/3,eta))
            + 3*libnest.units.T2/40*((2+libnest.units.T2X2/libnest.units.T2)
                                     *F_x(5/3,eta)
                                     -(1/2+libnest.units.T2X2/libnest.units.T2)
                                     *F_x(8/3,eta))*rho*kF**2
            + libnest.units.T3/48*np.power(rho,libnest.units.ALPHA+1)*
                (3-(1+2*libnest.units.X3)*eta**2)
            + 3*libnest.units.T4/40*kF**2*np.power(rho, libnest.units.BETA+1)
                *((2+libnest.units.X4)*F_x(5/3,eta)-(1/2+libnest.units.X4)*F_x(8/3,eta))
            + 3*libnest.units.T5/40*kF*np.power(rho,libnest.units.GAMMA+1)
                *((2+libnest.units.X5)*F_x(5/3,eta)-(1/2+libnest.units.X5)*F_x(8/3,eta)))



if __name__ == '__main__':
    pass

