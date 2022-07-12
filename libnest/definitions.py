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
    return libnest.units.hbar22M0/B*vsf + libnest.units.VUNIT*A/libnest.units.HBARC

def v_NV(B,j,rho,A):
    """Returns the velocity (mass velocity). It is adjusted to the entrainment
effects (definition by Nicolas Chamel Valentin Allard)."""
    return libnest.units.VUNIT*(libnest.units.hbar22M0/B*j/rho + A/libnest.units.HBARC)

def energy_per_nucleon(rho_n, rho_p):
    """Returns the energy per nucleon on infinite nuclear matter of given
    density of protons and neutrons, rho_p and rho_n, respectively, in MeV"""  
    kF = (3.**np.pi*(rho_n+rho_p))**(1./3.)
    eta = (rho_n-rho_p)/(rho_n+rho_p)
    rho = rho_n+rho_p
    F_x_5 = 0.5*((1+(rho_n-rho_p)/rho)**(5./3.)+(1-(rho_n-rho_p)/rho)**(5./3.))
    F_x_8 = 0.5*((1+(rho_n-rho_p)/rho)**(8./3.)+(1-(rho_n-rho_p)/rho)**(8./3.))    
    return (3*libnest.units.HBARC**2./20.*kF**2*((np.power(1+eta,5./3.)/
                                                   libnest.units.MN
                                                   +np.power(1-eta,5./3.)
                                                   /libnest.units.MP))
            + libnest.units.T0/8.*rho*(3-(2*libnest.units.X0 + 1)*eta**2) 
            + 3*libnest.units.T1/40*rho*kF**2*((2+libnest.units.X1)*F_x_5
                                               -(1./2.+libnest.units.X1)*F_x_8)
            + 3*libnest.units.T2/40*((2+libnest.units.T2X2/libnest.units.T2)
                                     *F_x_5-(1./2.+libnest.units.T2X2
                                             /libnest.units.T2)*F_x_8)*rho*kF**2
            + libnest.units.T3/48*np.power(rho,libnest.units.ALPHA+1)*
                (3-(1+2.*libnest.units.X3)*eta**2)
            + 3*libnest.units.T4/40.*kF**2*np.power(rho, libnest.units.BETA+1)
                *((2+libnest.units.X4)*F_x_5-(1./2.+libnest.units.X4)*F_x_8)
            + 3*libnest.units.T5/40*kF*np.power(rho,libnest.units.GAMMA+1)
                *((2+libnest.units.X5)*F_x_5-(1./2+libnest.units.X5)*F_x_8))

def effective_mass(rho, Ms, Mv):
    """Returns the effective mass of a nucleon of charge q given rho is the
    ratio of the density of matter of the specified type to total density.
    (rho = rho_n + rho_p)
    Assume Ms represents Ms/M, which is a unitless fraction, and Mv represents
    Mv/M, accordingly. ???
    For neutron matter, rho = rho_n
    For symmetric matter, rho = 2*rho_n"""
    return 1/(2*rho/Ms + (1-2*rho)/Mv)  

def q_effective_mass(M_q, rho, rho_q):
    """Returns effective mass Mq*/M of neutron or proton"""
    C_rho = libnest.units.T1/4*((1+libnest.units.X1/2)*rho-
                                (1./2+libnest.units.X1)*rho_q)
    +libnest.units.T4/4 * rho**libnest.units.BETA*((
        1+libnest.units.X4/2)*rho-(1./2*libnest.units.X4)*rho_q)
    +1./4*((libnest.units.T2+libnest.units.T2X2 /2)*rho 
          +(1./2*libnest.units.T2+libnest.units.T2X2)*rho_q)
    +libnest.units.T5/4*((1+libnest.units.X5/2)*rho+
                         (1./2+libnest.units.X5)*rho_q)*rho**libnest.units.GAMMA
    
    return (1./2*libnest.units.HBARC**2)/(((libnest.units.HBARC**2)/(2*M_q))
                                           +C_rho)

def symmetric_pairing_field(rho_n, rho_p):
    """Returns the pairing field for symmetric nuclear matter, for kF lower 
    than 1.38 fm^-1"""
    kF = (3.**np.pi*(rho_n+rho_p))**(1./3.)
    return 3.37968*(kF**2)*(kF-1.38236)**2/(kF**2+0.556092**2)/((kF-1.38236)**2
                                                                + 0.327517**2)

def neutron_pairing_field(rho_n):
    """Returns the pairing field for pure neutron matter, with kF lower than 
    1.31 fm^-1"""
    kF = kF = (3.**np.pi*(rho_n))**(1./3.)
    return 11.5586*(kF**2)*(kF-1.3142)**2/(kF**2+0.489932**2)/((kF-1.3142)**2
                                                                + 0.906146**2)

def neutron_ref_pairing_field(rho_n, rho_p):
    """Returns the reference pairing field for neutrons in uniform matter"""
    return (symmetric_pairing_field(rho_n, rho_p)*(1-abs((rho_n-rho_p)/
                                                        (rho_n+rho_p)))
            +neutron_pairing_field(rho_n)*rho_n/(rho_n+rho_p)*
            (rho_n-rho_p)/(rho_n+rho_p))

def proton_ref_pairing_field(rho_n, rho_p):
    """Returns the reference pairing field for protons in uniform matter"""
    return (symmetric_pairing_field(rho_n, rho_p)*(1-abs((rho_n-rho_p)/
                                                        (rho_n+rho_p)))
            -neutron_pairing_field(rho_n)*rho_n/(rho_n+rho_p)*
            (rho_n-rho_p)/(rho_n+rho_p))

def mean_field_potential(Mq, rho_n, rho_p, rho_q):
    """Returns the mean field potential
    rho_q is either rho_n or rho_p"""
    return libnest.units.HBARC**2/(2*Mq)+libnest.units.T1/4*(
        (1+libnest.units.X1/2)*(rho_n+rho_p)-(1./2+libnest.units.X1)*rho_q)
    +libnest.units.T4/4 * (rho_n+rho_p)**libnest.units.BETA*((
        1+libnest.units.X4/2)*(rho_n+rho_p)-(1./2*libnest.units.X4)*rho_q)
    +1./4*((libnest.units.T2+libnest.units.T2X2 /2)*(rho_n+rho_p) 
          +(1./2*libnest.units.T2+libnest.units.T2X2)*rho_q)
    +((1+libnest.units.X5/2)*(rho_n+rho_p)+(1./2+libnest.units.X5)*rho_q)
    *(rho_n+rho_p)**libnest.units.GAMMA*libnest.units.T5/4


#def isovector_effective_mass(M_q, rho, rho_q):
#    return


if __name__ == '__main__':
    pass
