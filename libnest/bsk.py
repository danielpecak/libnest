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
#================================
# TODO list:
# Formulas from https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804:
# Code formulas: A13, and dependent (A14, A15)
#

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
    return (3*units.HBARC**2/20*kF**2*((np.power(1+eta,5/3)/units.MN
                                                +np.power(1-eta,5/3)
                                                          /units.MP))
            + T0/8*rho*(3-(2*X0 + 1)*eta**2)
            + 3*T1/40*rho*kF**2*((2+X1)*F_x(5./3,eta)
                                               -(1/2+X1)*F_x(8./3,eta))
            + 3/40*((2*T2+T2X2)*F_x(5./3,eta)
                                     -(1./2*T2+T2X2)*F_x(8./3,eta))*rho*kF**2
            + T3/48*np.power(rho,ALPHA+1)*(3-(1+2*X3)*eta**2)
            + 3*T4/40*kF**2*np.power(rho, BETA+1)
                *((2+X4)*F_x(5./3,eta)-(1/2+X4)*F_x(8./3,eta))
            + 3*T5/40*kF*np.power(rho,GAMMA+1)
                *((2+X5)*F_x(5./3,eta)-(1/2+X5)*F_x(8./3,eta)))

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
