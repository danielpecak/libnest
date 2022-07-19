#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Description
# Original paper with formulas
# https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804
# NOTE: Table I is outdated; use the following data:
# forBSk31
# According to Goriely, Chamel, Pearson PRC 93 034337 (2016)
#================================
import sys
import numpy as np
from libnest import units
from libnest.units import HBARC, DENSEPSILON, NUMZERO
from libnest.units import MN, MP, HBAR2M_n, HBAR2M_p

T0   =-2302.01     # skyrme parameter t0 [MeV*fm^3]
T1   =762.99       # skyrme parameter t1 [MeV*fm<sup>5</sup>]
T2   =0.0          # skyrme parameter t2 [MeV*fm<sup>5</sup>]
T3   =13797.83     # skyrme parameter t3 [MeV*fm^(3+3*ALPHA)]
T4   =-500.         # skyrme parameter t4 [MeV*fm^(5+3*BETA)]
T5   =-40.          # skyrme parameter t5 [MeV*fm^(5+3*GAMMA)]
X0   =0.676655     # skyrme parameter x0 [1]
X1   =2.658109     # skyrme parameter x1 [1]
T2X2 =-422.29      # skyrme parameter x2t2 [1][MeV*fm<sup>5</sup>]
X3   =0.83982      # skyrme parameter x3 [1]
X4   =5.           # skyrme parameter x4 [1]
X5   =-12.         # skyrme parameter x5 [1]
ALPHA =(1./5.)     # [1]
BETA  =(1./12.)     # [1]
GAMMA =(1./4.)     # [1]
# NOTE: these are not important at the moment
YW   =2.           # [1]
FNP  =1.00         # [1]
FNM  =1.06         # [1]
FPP  =1.00         # [1]
FPM  =1.04         # [1]
KAPPAN =-36630.4   # [MeV*fm<sup>8</sup>]
KAPPAP =-45207.2   # [MeV*fm<sup>8</sup>]


# ================================
#       Auxiliary
# ================================
def rho2kf(rho):
    """Returns wavevector kF based on density rho."""
    return (3.*np.pi*np.pi*rho)**(1./3.)

def rho2tau(rho):
    """Returns kinetic density tau for uniform Fermi system of density rho."""
    return 0.6*(3.*np.pi)**(2./3.)*rho**(5./3.)

def rhoEta(rho_n, rho_p):
    """Returns total density and difference of densities from neutron and
    proton densities."""
    return rho_n+rho_p, rho_n-rho_p

# ================================
#       Pairing fields
# ================================
def symmetric_pairing_field(rho_n, rho_p):
    """Returns the pairing field for symmetric nuclear matter, for kF lower
    than 1.38 fm^-1.
    Formula (5.12) from NeST.pdf"""
    # NOTE: because of the IF statement to use this function for plotting eg.
    #   r = np.arange(0,0.2,0.0001)
    #   f = symmetric_pairing_field
    #   f=np.vectorize(f) # ADD THIS LINE TO MAKE IT WORK WITH NUMPY
    #   ax.plot(r,f(r),label='NeuM')
    kF = rho2kf((rho_n+rho_p)/2)
    delta = 3.37968*(kF**2)*((kF-1.38236)**2)/(((kF**2)+(0.556092**2))*
                                              ((kF-1.38236)**2+(0.327517**2)))
    i = np.where(kF>1.38)
    delta[i] = NUMZERO
    # if(kF > 1.38):
    #     return NUMZERO
    return delta

def neutron_pairing_field(rho_n):
    """Returns the pairing field for pure neutron matter, with kF lower than
    1.31 fm^-1.
    Formula (5.11) from NeST.pdf"""
    # NOTE: because of the IF statement to use this function for plotting eg.
    #   r = np.arange(0,0.2,0.0001)
    #   f = symmetric_pairing_field
    #   f=np.vectorize(f) # ADD THIS LINE TO MAKE IT WORK WITH NUMPY
    #   ax.plot(r,f(r),label='NeuM')
    kF = rho2kf(rho_n)
    delta = 11.5586*(kF**2)*((kF-1.3142)**2)/(((kF**2)+(0.489932**2))*
                                             (((kF-1.3142)**2)+(0.906146**2)))
    i = np.where(kF>1.31)
    delta[i] = NUMZERO
    # if(kF > 1.31):
    #     return NUMZERO
    return delta

def neutron_ref_pairing_field(rho_n, rho_p):
    """Returns the reference pairing field for neutrons in uniform matter.
    Formula (5.10) from NeST.pdf"""
    rho, eta = rhoEta(rho_n, rho_p)
    rho = rho + DENSEPSILON
    return (symmetric_pairing_field(rho_n, rho_p)*(1-abs(eta/rho))
            +neutron_pairing_field(rho_n)*rho_n/rho*eta/rho)

def proton_ref_pairing_field(rho_n, rho_p):
    """Returns the reference pairing field for protons in uniform matter.
    Formula (5.10) from NeST.pdf"""
    # TODO: use these definitions
    rho, eta = rhoEta(rho_n, rho_p) #eta = rho_n - rho_p
    rho = rho + DENSEPSILON
    return (symmetric_pairing_field(rho_n, rho_p)*(1-abs(eta/rho))
            -neutron_pairing_field(rho_n)*rho_p/rho*eta/rho)



# ================================
#       Effective masses
# ================================
# TODO list:
# https://journals.aps.org/prc/pdf/10.1103/PhysRevC.82.035804
# Code function: Eq. 10 and plot it as a function of rho (I suppose this is what is done in
# Fig 5 in the inset). q=n or p (neutron or proton counterpart)
# rho = rho_n + rho+p
# Make plots for neutron matter (NeuM) where rho = rho_n [rho_p =0 ]
# Make plots for symmetric matter (SM) where rho = 2*rho_n [rho_p =rho_n ]
#
def effective_mass(rho, Ms, Mv):
    """Returns the effective mass of a nucleon of charge q given rho is the
    ratio of the density of matter of the specified type to total density.
    (rho = rho_n + rho_p)
    Assume Ms represents Ms/M, which is a unitless fraction, and Mv represents
    Mv/M, accordingly. ???
    For neutron matter, rho = rho_n
    For symmetric matter, rho = 2*rho_n"""
    # NOTE: I think this probably would not be needed anymore
    return 1/(2*rho/Ms + (1-2*rho)/Mv)

def isoscalarM(rho_n, rho_p):
    """Calculated effective isoscalar mass M_s for a given uniform system
    with neutron and proton densities rho_n, rho_p respectively.  """
    # TODO FILL IN using effMn and effMp
    pass

def isovectorM(rho_n, rho_p):
    """Calculated effective isovector mass M_s for a given uniform system
    with neutron and proton densities rho_n, rho_p respectively.  """
    # TODO FILL IN using effMn and effMp
    pass

def effMn(rho_n, rho_p):
    """Effective mass of a neutron in nuclear medium."""
    # TODO please test it: make plots as a function of density x
    return HBAR2M_n/B_q(rho_n, rho_p,'n')

def effMp(rho_n, rho_p):
    """Effective mass of a proton in nuclear medium."""
    # TODO please test it: make plots as a function of density x
    return HBAR2M_p/B_q(rho_n, rho_p,'p')

def U_q(rho_n, rho_p,q): # TODO
    """Returns the mean field potential B_q
    rho_q is either rho_n or rho_p.
    Formula (5.14) from NeST.pdf
    """
    # TODO: modify IF statement if needed
    if(q=='n'):
        rho_q = rho_n
        rho_q_prime = rho_p
        HBAR2M_q = HBAR2M_n
    elif(q=='p'):
        rho_q = rho_p
        rho_q_prime = rho_n
        HBAR2M_q = HBAR2M_p
    else:
        sys.exit('# ERROR: Nucleon q must be either n or p')
    rho = rho_n + rho_p
    # TODO fill in what is returned
    return T0*((1+0.5*X0)*rho-(0.5+X0)*rho_q)+T3/12.*np.power(rho,(ALPHA-1))*(
        ((0.5+0.5*X3)*rho**2*(ALPHA+2))-(0.5+X3)*(2*rho*rho_q*ALPHA*
                                                  (rho_q_prime)**2))


def B_q(rho_n, rho_p, q):
    """Returns the mean field potential B_q
    rho_q is either rho_n or rho_p.
    Formula (5.13) from NeST.pdf
    """
    if(q=='n'):
        rho_q = rho_n
        HBAR2M_q = HBAR2M_n
    elif(q=='p'):
        rho_q = rho_p
        HBAR2M_q = HBAR2M_p
    else:
        sys.exit('# ERROR: Nucleon q must be either n or p')
    rho = rho_n + rho_p
    return (HBAR2M_q
             + T1/4.*((1.+X1/2.)*rho - (1./2.+X1)*rho_q)
             + T4/4.*((1.+X4/2.)*rho - (1./2.+X4)*rho_q)*np.power(rho,BETA)
             + 1./4.*((T2+T2X2/2.)*rho+(1./2.*T2+T2X2)*rho_q)
             + T5/4.*((1.+X5/2.)*rho + (1./2.+X5)*rho_q)*np.power(rho,GAMMA)
             )

# TODO list:
# Formulas from https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804:
# Code formulas: A13, and dependent (A14, A15)
#



def energy_per_nucleon(rho_n, rho_p):
    """Returns the energy per nucleon on infinite nuclear matter of given
    density of protons and neutrons, rho_p and rho_n, respectively, in MeV.
    Formula (A13) from https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804
    """
    # TODO: testing: making plots for total density [0, 0.2] and [0, 1.0]
    #               for pure neutron and symmetric nuclear matter X
    rho = rho_n+rho_p+DENSEPSILON
    kF = rho2kf(0.5*rho) # Formula (A14) from the paper is using rho/2
    eta = (rho_n-rho_p)/rho
    F_x_5 = 0.5*((1+eta)**(5./3.)+(1-eta)**(5./3.)) # Formula (A15)
    F_x_8 = 0.5*((1+eta)**(8./3.)+(1-eta)**(8./3.)) # Formula (A15)

    return (3*(HBARC**2)/20*(kF**2)*((np.power((1+eta),5/3)/MN
                                                +np.power((1-eta),5/3)/MP))
            + T0/8.*rho*(3-(2*X0 + 1)*eta**2)
            + 3.*T1/40*rho*(kF**2)*((2+X1)*F_x_5 -(0.5+X1)*F_x_8)
            + 3./40*((2*T2+T2X2)*F_x_5 +(1./2*T2+T2X2)*F_x_8)*rho*(kF**2)
            + T3/48*np.power(rho,(ALPHA+1))*(3-(1+2*X3)*(eta**2))
            + 3*T4/40*(kF**2)*np.power(rho,BETA+1)*((2+X4)*F_x_5-(1/2+X4)*F_x_8)
            + 3*T5/40*(kF**2)*np.power(rho,GAMMA+1)*((2+X5)*F_x_5+(1/2+X5)*F_x_8))



# TODO list:
# PHYSICAL REVIEW C 104, 055801 (2021)
# NOTE: densities such as TAU, MU, J will be given in the future from the data
# Formulas 9-14,23, A8-A10
