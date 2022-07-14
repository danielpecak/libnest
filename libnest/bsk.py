#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Description
# Original paper with formulas
# https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804
# NOTE: Table I is outdated; use the following data:
# forBSk31
# According to Goriely, Chamel, Pearson PRC 93 034337 (2016)
#================================
from libnest import units
import numpy as np

T0   =-2302.01     # skyrme parameter t0 [MeV*fm^3]
T1   =762.99       # skyrme parameter t1 [MeV*fm<sup>5</sup>]
T2   =0.0          # skyrme parameter t2 [MeV*fm<sup>5</sup>]
T3   =13797.83     # skyrme parameter t3 [MeV*fm^(3+3*ALPHA)]
T4   =-500         # skyrme parameter t4 [MeV*fm^(5+3*BETA)]
T5   =-40          # skyrme parameter t5 [MeV*fm^(5+3*GAMMA)]
X0   =0.676655     # skyrme parameter x0 [1]
X1   =2.658109     # skyrme parameter x1 [1]
T2X2 =-422.29      # skyrme parameter x2t2 [1][MeV*fm<sup>5</sup>]
X3   =0.83982      # skyrme parameter x3 [1]
X4   =5.           # skyrme parameter x4 [1]
X5   =-12.         # skyrme parameter x5 [1]
ALPHA =(1./5.)     # [1]
BETA =(1./12.)     # [1]
GAMMA =(1./4.)     # [1]
#define YW  2.       // [1]
#define FNP 1.00     // [1]
#define FNM 1.06     // [1]
#define FPP 1.00     // [1]
#define FPM 1.04     // [1]
#define KAPPAN -36630.4 // [MeV*fm<sup>8</sup>]
#define KAPPAP -45207.2 // [MeV*fm<sup>8</sup>]

# ================================
#       Auxiliary
# ================================
def rho2kf(rho):
    """Returns wavevector kF based on density rho."""
    return (3.**np.pi*rho)**(1./3.)

# ================================
#       Pairing fields
# ================================
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

# ================================
#       Effective masses
# ================================
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
    C_rho = units.T1/4*((1+units.X1/2)*rho-
                                (1./2+units.X1)*rho_q)
    +units.T4/4 * rho**units.BETA*((
        1+units.X4/2)*rho-(1./2*units.X4)*rho_q)
    +1./4*((units.T2+units.T2X2 /2)*rho
          +(1./2*units.T2+units.T2X2)*rho_q)
    +units.T5/4*((1+units.X5/2)*rho+
                         (1./2+units.X5)*rho_q)*rho**units.GAMMA

    return (1./2*units.HBARC**2)/(((units.HBARC**2)/(2*M_q))
                                           +C_rho)

def mean_field_potential(Mq, rho_n, rho_p, rho_q):
    """Returns the mean field potential
    rho_q is either rho_n or rho_p"""
    return units.HBARC**2/(2*Mq)+units.T1/4*(
        (1+units.X1/2)*(rho_n+rho_p)-(1./2+units.X1)*rho_q)
    +units.T4/4 * (rho_n+rho_p)**units.BETA*((
        1+units.X4/2)*(rho_n+rho_p)-(1./2*units.X4)*rho_q)
    +1./4*((units.T2+units.T2X2 /2)*(rho_n+rho_p)
          +(1./2*units.T2+units.T2X2)*rho_q)
    +(((1+units.X5/2)*(rho_n+rho_p)+(1./2+units.X5)*rho_q)*
    (rho_n+rho_p)**units.GAMMA*units.T5/4)

# TODO list:
# Formulas from https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804:
# Code formulas: A13, and dependent (A14, A15)
#



def energy_per_nucleon(rho_n, rho_p):
    """Returns the energy per nucleon on infinite nuclear matter of given
    density of protons and neutrons, rho_p and rho_n, respectively, in MeV"""
    # TODO: kF and rho are connected see kF(rho)
    # TODO: F_x not needed as an argument
    # TODO: rho_n and rho_p as an argument, instead of rho and eta
    # NOTE: removed '/T2' since it's division by 0
    # NOTE: 2/3 in many programming languages gives 0 (operation on integeres)
    #       therefore I prefer to write explicitely 2./3 that would work in
    #       C or FORTRAN - just in case
    rho = rho_n+rho_p+units.DENSEPSILON
    kF = (3.*np.pi**2*rho/2.)**(1./3.)
    eta = (rho_n-rho_p)/rho
    F_x_5 = 0.5*((1+eta)**(5./3.)+(1-eta)**(5./3.))
    F_x_8 = 0.5*((1+eta)**(8./3.)+(1-eta)**(8./3.))

    return (3*units.HBARC**2/20*kF**2*((np.power(1+eta,5/3)/units.MN
                                                +np.power(1-eta,5/3)
                                                          /units.MP))
            + T0/8.*rho*(3-(2*X0 + 1)*eta**2)
            + 3.*T1/40*rho*kF**2*((2+X1)*F_x_5 -(1/2+X1)*F_x_8)
            + 3./40*((2*T2+T2X2)*F_x_5 +(1./2*T2+T2X2)*F_x_8)*rho*kF**2
            + T3/48*np.power(rho,ALPHA+1)*(3-(1+2*X3)*eta**2)
            + 3*T4/40*kF**2*np.power(rho, BETA+1)*((2+X4)*F_x_5-(1/2+X4)*F_x_8)
            + 3*T5/40*kF*np.power(rho,GAMMA+1)*((2+X5)*F_x_5+(1/2+X5)*F_x_8))

# TODO (5.13) from NeST.pdf

# TODO (5.11) from NeST.pdf
# TODO (5.12) from NeST.pdf
# TODO (5.10) from NeST.pdf

# TODO list:
# https://journals.aps.org/prc/pdf/10.1103/PhysRevC.82.035804
# Code function: Eq. 10 and plot it as a function of rho (I suppose this is what is done in
# Fig 5 in the inset). q=n or p (neutron or proton counterpart)
# rho = rho_n + rho+p
# Make plots for neutron matter (NeuM) where rho = rho_n [rho_p =0 ]
# Make plots for symmetric matter (SM) where rho = 2*rho_n [rho_p =rho_n ]
#
# TODO list:
# PHYSICAL REVIEW C 104, 055801 (2021)
# NOTE: densities such as TAU, MU, J will be given in the future from the data
# Formulas 9-14,23, A8-A10
