#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Description
# Original paper with formulas
# https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804
# NOTE: Table I is outdated; use the following data:
# forBSk31
# According to Goriely, Chamel, Pearson PRC 93 034337 (2016)
#================================
"""
BSk.py
======
The module that contains Brussels-Montreal parametrization with formulas
for energy density, effective masses, etc.
Both for uniform system and general expressions.


The constants for BSk31:

.. data:: T0   =-2302.01

    Skyrme parameter :math:`t_0` [MeV fm :sup:`3`]

.. data:: T1   =762.99

    Skyrme parameter :math:`t_1` [MeV fm :sup:`5`]

.. data:: T2   =0.0

    Skyrme parameter :math:`t_2` [MeV fm :sup:`5`]

.. data:: T3   =13797.83

    Skyrme parameter :math:`t_3` [MeV fm :sup:`(3+3*ALPHA)`]

.. data:: T4   =-500.

    Skyrme parameter :math:`t_4` [MeV fm :sup:`(5+3*BETA)`]

.. data:: T5   =-40.

    Skyrme parameter :math:`t_5` [MeV fm :sup:`(5+3*GAMMA)`]

.. data:: X0   =0.676655

    Skyrme parameter :math:`x_0` [1]

.. data:: X1   =2.658109

    Skyrme parameter :math:`x_1` [1]

.. data:: T2X2 =-422.29

    Skyrme parameter :math:`x_2t_2` [1][MeV fm :sup:`5`]

.. data:: X3   =0.83982

    Skyrme parameter :math:`x_3` [1]

.. data:: X4   =5.

    Skyrme parameter :math:`x_4` [1]

.. data:: X5   =-12.

    Skyrme parameter :math:`x_5` [1]

.. data:: ALPHA =(1./5.)

    [1]

.. data:: BETA  =(1./12.)

    [1]

.. data:: GAMMA =(1./4.)

    [1]

Note:
    these are not important at the moment

.. data:: YW   =2.

    [1]

.. data:: FNP  =1.00

    [1]

.. data:: FNM  =1.06

    [1]

.. data:: FPP  =1.00

    [1]

.. data:: FPM  =1.04

    [1]

.. data:: KAPPAN =-36630.4

    [MeV fm :sup:`8`]

.. data:: KAPPAP =-45207.2

    [MeV fm :sup:`8`]


==========

"""
import sys
import numpy as np
from libnest import units
from libnest.units import HBARC, DENSEPSILON, NUMZERO
from libnest.units import MN, MP, HBAR2M_n, HBAR2M_p
from libnest.definitions import rho2kf, rhoEta

T0   =-2302.01     # Skyrme parameter :math:`t_0` [MeV fm :sup:`3`]
T1   =762.99       # Skyrme parameter :math:`t_1` [MeV*fm :sup:`5`]
T2   =0.0          # skyrme parameter t2 [MeV*fm<sup>5</sup>]
T3   =13797.83     # skyrme parameter t3 [MeV*fm^(3+3*ALPHA)]
T4   =-500.        # skyrme parameter t4 [MeV*fm^(5+3*BETA)]
T5   =-40.         # skyrme parameter t5 [MeV*fm^(5+3*GAMMA)]
X0   =0.676655     # skyrme parameter x0 [1]
X1   =2.658109     # skyrme parameter x1 [1]
T2X2 =-422.29      # skyrme parameter x2t2 [1][MeV*fm<sup>5</sup>]
X3   =0.83982      # skyrme parameter x3 [1]
X4   =5.           # skyrme parameter x4 [1]
X5   =-12.         # skyrme parameter x5 [1]
ALPHA =(1./5.)     # [1]
BETA  =(1./12.)    # [1]
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
#       Pairing fields
# ================================
def neutron_pairing_field(rho_n):
    #    Formula (5.12) from NeST.pdf
    """
    Returns the pairing field for neutron nuclear matter for kF lower
    than 1.38 fm^-1 (where the field is zero).
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

    Returns:
        float: pairing field for symmetric matter :math:`\\Delta_{sym}` [fm :sup:`-3`]
    
    """
    kF = rho2kf((rho_n)/2)
    delta = 3.37968*(kF**2)*((kF-1.38236)**2)/(((kF**2)+(0.556092**2))*
                                              ((kF-1.38236)**2+(0.327517**2)))
    i = np.where(kF>1.38)
    delta[i] = NUMZERO
    #delta = np.ma.masked_invalid(delta)
    return delta

def symmetric_pairing_field(rho_n, rho_p):
    #   Formula (5.11) from NeST.pdf
    """
    Returns the pairing field for pure symmetric matter for kF lower than
    1.31 fm^-1 (where the field is zero).
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        
    Returns:
        float: pairing field for neutron matter :math:`\\Delta_{NeuM}` [fm :sup:`-3`]
    """
    kF = rho2kf((rho_n+rho_p)/2)
    delta = 11.5586*(kF**2)*((kF-1.3142)**2)/(((kF**2)+(0.489932**2))*
                                             (((kF-1.3142)**2)+(0.906146**2)))
    i = np.where(kF>1.31)
    delta[i] = NUMZERO

    return delta

def neutron_ref_pairing_field(rho_n, rho_p):
    #   Formula (5.10) from NeST.pdf
    """
    Returns the reference pairing field for neutrons in uniform matter.
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        
    Returns:
        float: pairing field for neutrons :math:`\\Delta_n` [fm :sup:`-3`]
    """
    rho, eta = rhoEta(rho_n, rho_p)
    rho = rho + DENSEPSILON
    return (symmetric_pairing_field(rho_n, rho_p)*(1-abs(eta/rho))
            +neutron_pairing_field(rho_n)*rho_n/rho*eta/rho)

def proton_ref_pairing_field(rho_n, rho_p):
    #   Formula (5.10) from NeST.pdf
    """
    Returns the reference pairing field for protons in uniform matter.
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        
    Returns:
        float: pairing field for protons :math:`\\Delta_p` [fm :sup:`-3`]
    """
    rho, eta = rhoEta(rho_n, rho_p) #eta = rho_n - rho_p
    rho = rho + DENSEPSILON
    return (symmetric_pairing_field(rho_n, rho_p)*(1-abs(eta/rho))
            -neutron_pairing_field(rho_n)*rho_p/rho*eta/rho)

def eF_n(kF):
    """
    Returns Fermi energy for neutrons based on wavevector kF.

    Args:
        kF (float):  wavevector :math:`k_F`

    Returns:
        float: Fermi energy :math:`\\Epsilon_F` [MeV]
    """
    return HBAR2M_n * kF**2

def E_minigap_n(rho_n):
    """
    Returns the energy of minigap :math:`E_{mg}` [MeV] for neutron matter.
    
    Parameters
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

    Returns
        float: energy of minigap :math:`E_{mg}` [MeV]
    """
    delta = neutron_ref_pairing_field(rho_n, 0.)
    return 4./3. * np.abs(delta)**2/eF_n(rho2kf(rho_n))


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

def isoscalarM(rho_n, rho_p):
    """
    Calculates effective isoscalar mass M_s for a given uniform system
    with neutron and proton densities rho_n and rho_p respectively. 
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        
    Returns:
        float: effective isoscalar mass :math:`M_{s}^{*}` [MeV]
        
    See also:
        :func:`effMn`
        :func:`effMp`
    
    """
    return 2/((1./effMn(rho_n, rho_p))+(1./effMp(rho_n, rho_p)))

def isovectorM(rho_n, rho_p):
    """
    Calculated effective isovector mass :math:`M_v` for a given uniform system
    with neutron and proton densities rho_n, rho_p respectively.
    
     Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        
    Returns:
        float: effective isovector mass :math:`M_{v}^{*}` [MeV]
    
    See also:
        :func:`effMn`
        :func:`effMp`
    """
    Mp = effMp(rho_n, rho_p)
    rho = rho_n + rho_p + DENSEPSILON
    rho = rho + DENSEPSILON
    return isoscalarM(rho_n, rho_p)*Mp*(2*rho_p-rho)/(2*Mp*rho_p-rho)

def effMn(rho_n, rho_p):
    """
    Returns the effective mass of a neutron in nuclear medium.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

    Returns:
        float: effective mass of a neutron :math:`M_{n}^{*}` [MeV]

    See also:
        :func:`B_q`
    """
    return HBAR2M_n/B_q(rho_n, rho_p,'n')

def effMp(rho_n, rho_p):
    """
    Returns the effective mass of a proton in nuclear medium.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

    Returns:
        float: effective mass of a proton :math:`M_{p}^{*}` [MeV]

    See also:
        :func:`B_q`
    """
    return HBAR2M_p/B_q(rho_n, rho_p,'p')

def U_q(rho_n, rho_p,q):
    #     Formula (5.14) from NeST.pdf
    """
    Returns the mean field potential from density :math:`\\rho` variation.

    rho_q is either rho_n or rho_p.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        q (string): nucleon type choice ('p' - proton, or 'n' - neutron)

    Returns:
        float: Mean field potential :math:`U_q` [MeV]
    """
    if(q=='n'):
        rho_q = rho_n + DENSEPSILON
        rho_q_prime = rho_p + DENSEPSILON
    elif(q=='p'):
        rho_q = rho_p + DENSEPSILON
        rho_q_prime = rho_n + DENSEPSILON
    else:
        sys.exit('# ERROR: Nucleon q must be either n or p')
    rho = rho_n + rho_p + DENSEPSILON
    return T0*((1+0.5*X0)*rho-(0.5+X0)*rho_q)+T3/12.*np.power(rho,(ALPHA-1))*(
        ((0.5+0.5*X3)*rho**2*(ALPHA+2))-(0.5+X3)*(2*rho*rho_q*ALPHA*
                                                  (rho_q_prime)**2))

def B_q(rho_n, rho_p, q):
    #    Formula (5.13) from NeST.pdf
    """
    Returns the mean field potential B_q (coming from variation over kinetic
    density, or effective mass). rho_q is either rho_n or rho_p.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        q (string): nucleon type choice ('p' - proton, or 'n' - neutron)

    Returns:
        float: effective mass of a proton :math:`M_{p}^{*}` [MeV fm :sup:`2`]
    """
    if(q=='n'):
        rho_q = rho_n + DENSEPSILON
        HBAR2M_q = HBAR2M_n
    elif(q=='p'):
        rho_q = rho_p + DENSEPSILON
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

# Formulas from https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804:
# Code formulas: A13, and dependent (A14, A15)
def energy_per_nucleon(rho_n, rho_p):
    """
    Returns the energy per nucleon on infinite nuclear matter of given
    density of protons and neutrons, rho_p and rho_n, respectively, in MeV.
    Formula (A13) from https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804
    :cite:p:`chamel2009further`.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

    Returns:
        float: energy per neutron :math:`e_n` [MeV]
    """
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


# ================================
#        Energy functional
# ================================
# TODO list:
# PHYSICAL REVIEW C 104, 055801 (2021)
# NOTE: densities such as TAU, MU, J will be given in the future from the data
# Formulas 9-14,23, A8-A10

def C_rho(rho):
    """
    Calculates the energy functional :math:`C^{\\rho}` cooefficient for NeuM

    Args:
        rho (float): particle density :math:`\\rho` [fm :sup:`-3`]

    Returns:
        float: C coefficient :math:`C^{\\rho}`
    """
    return -1./4.*T0*(X0-1)-T3/24.*(X3-1)*rho**ALPHA

def C_tau(rho):
    """
    Calculates the energy functional :math:`C^{\\tau}` cooefficient for NeuM

    Args:
        rho (float): particle density :math:`\\rho` [fm :sup:`-3`]

    Returns:
        float: C coefficient :math:`C^{\\tau}`
    """
    return (-T1/8.*(X1-1.)+3./8.*T2X2 +3./8.*T2-T4/8.*(X4-1.)*rho**BETA+3./8.*
            T5*(X5+1.)*rho**GAMMA)

def C_delta_rho(rho):
    """
    Calculates the energy functional :math:`C^{\\Delta * \\rho}` cooefficient
    for NeuM

    Args:
        rho (float): particle density :math:`\\rho` [fm :sup:`-3`]

    Returns:
        float: C coefficient :math:`C^{\\Delta * \\rho}`
    """
    return (3./32.*T1*(X1-1.)+(3./32.*T2X2+3./32.*T2)+3./32.*T4*(X4-1.)*rho**BETA
             +3./32.*T5*(X5+1.)*rho**GAMMA)

def C0_rho(rho_n, rho_p):
    """
    Calculates the energy functional :math:`C^{\\rho}_0` cooefficient

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

    Returns:
        float: :math:`C^{\\rho}_0` coefficient
    """
    return 3. / 8. * T0 + 3. / 48. * T3 * np.power(rho_n + rho_p, ALPHA)


def C1_rho(rho_n, rho_p):
    """
    Calculates the energy functional :math:`C^{\\rho}_1` cooefficient

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

    Returns:
        float: :math:`C^{\\rho}_1` coefficient
    """
    return -1. / 4. * T0 * (1. / 2. + X0) - 1. / 24. * T3 * (1. / 2. + X3) * np.power(rho_n + rho_p, ALPHA)


def C0_tau(rho_n, rho_p):
    """
    Calculates the energy functional :math:`C^{\\tau}_0` cooefficient

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

    Returns:
        float: :math:`C^{\\tau}_0` coefficient
    """
    rho = rho_n + rho_p
    return (3. / 16.) * T1 + 5. / 16. * T2 + 1. / 4. * T2X2 + 3. / 16. * T4 * np.power(rho, BETA) + 1. / 4. * T5 * (1.25 + X5) * np.power(rho, GAMMA)


def C1_tau(rho_n, rho_p):
    """
    Calculates the energy functional :math:`C^{\\tau}_1` cooefficient

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

    Returns:
        float: :math:`C^{\\tau}_1` coefficient
    """
    rho = rho_n + rho_p
    return (-1. / 8.) * T1 * (1. / 2. + X1) + 1. / 16. * T2 + 1. / 8. * T2X2 - 1. / 8. * T4 * np.power(rho, BETA) * (1. / 2. + X4) + 1. / 8. * T5 * np.power(rho, GAMMA) * (1. / 2. + X5)

def epsilon_rho_np(rho_n, rho_p):
    """
    Energy functional :math:`\\epsilon_{\\rho}` for particle matter, related to
    the interaction of the nucleons with the background matter density.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

    Returns:
        float: energy functional :math:`\\epsilon_{\\rho}`
    """
    rho = rho_n + rho_p + DENSEPSILON
    return C0_rho(rho_n, rho_p) * np.power(rho, 2) + C1_rho(rho_n, rho_p) * np.power(rho_n - rho_p, 2)


def epsilon_tau_np(rho_n, rho_p, tau_n, tau_p, jsum2, jdiff2):
    """
    Energy functional :math:`\\epsilon_{\\tau}` for particle matter,
    related to the density-dependent effective mass. It gives rise to current couplings.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        tau_n (float): kinetic density :math:`\\tau` [fm :sup:`-5`]
        tau_n (float): kinetic density :math:`\\tau` [fm :sup:`-5`]
        jsum2 (float): sum of momentum density/current vectors :math:`j` [fm :sup:`-3`]
        jdiff2 (float) : difference of momentum density/current vectors :math:`j` [fm :sup:`-3`]

    Returns:
        float: energy functional :math:`\\epsilon_{\\tau}`
    """
    return C0_tau(rho_n, rho_p) * ((rho_n + rho_p) * (tau_n + tau_p) - jsum2) + C1_tau(rho_n, rho_p) * ((rho_n - rho_p) * (tau_n - tau_p) - jdiff2)


def epsilon_delta_rho_np(rho_n, rho_p, rho_grad_n_square, rho_grad_p_square, rho_grad_square):
    """
    Energy functional :math:`\\epsilon_{\\Delta \\rho}` for particle matter,
    related to the interaction of the nucleons with the background fluctuations.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        rho_grad_n_square (float): neutron density gradient squared :math:`\\nabla \\rho` [fm :sup:`-8`]
        rho_grad_p_square (float): proton density gradient squared :math:`\\nabla \\rho` [fm :sup:`-8`]
        rho_grad_square (float): total density gradient squared :math:`\\nabla \\rho` [fm :sup:`-8`]

    Returns:
        float: energy functional :math:`\\epsilon_{\\Delta \\rho}`
    """
    rho = rho_n+rho_p + DENSEPSILON
    grad_rho_n_rho_p = 0.5*(rho_grad_p_square-rho_grad_n_square-rho_grad_p_square)
  
    return (3./16.*T1*((1.+0.5*X1)*rho_grad_square-(0.5+X1)*(rho_grad_n_square+rho_grad_p_square))
            -1./16.*T2X2*(0.5*rho_grad_square+rho_grad_n_square+rho_grad_p_square)
            -1./16.*T2*(rho_grad_square+0.5*rho_grad_n_square+0.5*rho_grad_p_square)
            +3./16.*T4*np.power(rho,BETA)*((1.+0.5*X4)*rho_grad_square-(0.5+X4)*(rho_grad_n_square+rho_grad_p_square))
            -1./16.*T5*np.power(rho,GAMMA)*((1.+0.5*X5)*rho_grad_square+(0.5+X5)*(rho_grad_n_square+rho_grad_p_square))
            +BETA/8.*T4*np.power(rho,BETA-1.)*((1.+0.5*X4)*rho*rho_grad_square-(0.5+X4)
                                          *(rho_n*rho_grad_n_square+rho_p*rho_grad_p_square+rho*grad_rho_n_rho_p))
            )

def epsilon(rho_n, rho_p, rho_grad, tau, j, nu, q, kappa):
    """
    Calculates the total energy functional :math:`\\epsilon`, with a limit for
    single-particle energies :math:`\\epsilon_{\\Lambda}` = 6.5 MeV.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        rho_grad (float): particle density gradient :math:`\\nabla \\rho` [fm :sup:`-4`]
        tau (float): kinetic density :math:`\\tau` [fm :sup:`-5`]
        j (float): momentum density/current :math:`j` [fm :sup:`-3`]
        nu (float): anomalous density :math:`\\nu` [fm :sup:`-3`]
        q (string): nucleon type choice ('p' - proton, or 'n' - neutron)
        kappa (float):
            what is kappa? (no Eq.9 in Ref.41)

    Returns
        float: nergy functional :math:`\\epsilon`
    """
    rho = rho_n + rho_p + DENSEPSILON
    if(q=='n'):
        M = MN
    elif(q=='p'):
        M = MP
    else:
        sys.exit('# ERROR: Nucleon q must be either n or p')
    return HBARC**2/2./M*tau + epsilon_rho(rho)+epsilon_delta_rho(rho, rho_grad)+epsilon_tau(rho, tau, j)+epsilon_pi(rho_n, rho_p, rho_grad, nu, q, kappa)

def epsilon_np(rho_n, rho_p, rho_grad_n, rho_grad_p, tau_n, tau_p, jsum2, jdiff2, nu_n, nu_p, kappa_n, kappa_p):
    """
    Calculates the total energy functional :math:`\\epsilon`, with a limit for
    single-particle energies :math:`\\epsilon_{\\Lambda}` = 6.5 MeV.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        rho_grad_n (float): neutron density gradient :math:`\\nabla \\rho` [fm :sup:`-4`]
        rho_grad_p (float): proton density gradient :math:`\\nabla \\rho` [fm :sup:`-4`]
        tau_n (float): kinetic density :math:`\\tau` [fm :sup:`-5`]
        tau_n (float): kinetic density :math:`\\tau` [fm :sup:`-5`]
        jsum2 (float): sum of momentum density/current vectors :math:`j` [fm :sup:`-3`]
        jdiff2 (float) : difference of momentum density/current vectors :math:`j` [fm :sup:`-3`]
        nu_n (float): neutron anomalous density :math:`\\nu` [fm :sup:`-3`]
        nu_p (float): proton anomalous density :math:`\\nu` [fm :sup:`-3`]
        kappa_n (float):
        kappa_p (float):
            what is kappa? (no Eq.9 in Ref.41)
 
    TO DO: write eq for jsum2 and jdiff2

    Returns
        float: nergy functional :math:`\\epsilon`
    """
    rho_grad_n_square = (rho_grad_n)**2
    rho_grad_p_square = (rho_grad_p)**2
    rho_grad_square = rho_grad_n_square + rho_grad_p_square
    return (HBARC**2/2./MN*tau_n + HBARC**2/2./MP*tau_p
            + epsilon_rho_np(rho_n, rho_p)
            + epsilon_delta_rho_np(rho_n, rho_p, rho_grad_n_square, rho_grad_p_square, rho_grad_square)
            + epsilon_tau_np(rho_n, rho_p, tau_n, tau_p, jsum2, jdiff2)
            + epsilon_pi_np(rho_n, rho_p, rho_grad_n, rho_grad_p, nu_n, nu_p, kappa_n, kappa_p)
            )

def epsilon_rho(rho):
    """
    Energy functional :math:`\\epsilon_{\\rho}` for particle matter, related to
    the interaction of the nucleons with the background matter density.

    Args:
        rho (float): particle density :math:`\\rho` [fm :sup:`-3`]

    Returns:
        float: energy functional :math:`\\epsilon_{\\rho}`
    """
    return C_rho(rho)*rho**2
    
def epsilon_delta_rho(rho, rho_grad):
    """
    Energy functional :math:`\\epsilon_{\\Delta \\rho}` for particle matter,
    related to the interaction of the nucleons with the background fluctuations.

    Args:
        rho (float): particle density :math:`\\rho` [fm :sup:`-3`]
        rho_grad (float): particle density gradient :math:`\\nabla \\rho` [fm :sup:`-4`]

    Returns:
        float: energy functional :math:`\\epsilon_{\\Delta \\rho}`
    """
    return (-(rho_grad)**2)*C_delta_rho(rho)

def epsilon_tau(rho, tau, j):
    """
    Energy functional :math:`\\epsilon_{\\tau}` for particle matter,
    related to the density-dependent effective mass. It gives rise to current couplings.

    Args:
        rho (float): particle density :math:`\\rho` [fm :sup:`-3`]
        tau (float): kinetic density :math:`\\tau` [fm :sup:`-5`]
        j (float): momentum density/current :math:`j` [fm :sup:`-3`]

    Returns:
        float: energy functional :math:`\\epsilon_{\\tau}`
    """
    return (rho*tau-j**2)*C_tau(rho)

def epsilon_pi(rho_n, rho_p, rho_grad, nu, q, kappa):
    """
    Energy functional :math:`\\epsilon_{\\pi}` for particle matter,
    related to pairng energy density.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        rho_grad (float): particle density gradient :math:`\\nabla \\rho` [fm :sup:`-4`]
        nu (float): anomalous density :math:`\\nu` [fm :sup:`-3`]
        q (string): nucleon type choice ('p' - proton, or 'n' - neutron)
        kappa (float):
           - is this tern supposed to be included? 

    Returns:
        float: energy functional :math:`\\epsilon_{\\pi}`
    """
    return 1./4.*v_pi(rho_n, rho_p, q)*nu*(1-nu)+kappa*(np.abs(rho_grad))**2

def epsilon_pi_np(rho_n, rho_p, rho_grad_n, rho_grad_p, nu_n, nu_p, kappa_n, 
                    kappa_p):
    """
    Energy functional :math:`\\epsilon_{\\pi}` for particle matter,
    related to pairng energy density.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        rho_grad_n (float): neutron density gradient :math:`\\nabla \\rho` [fm :sup:`-4`]
        rho_grad_p (float): proton density gradient :math:`\\nabla \\rho` [fm :sup:`-4`]
        nu_n (float): neutron anomalous density :math:`\\nu` [fm :sup:`-3`]
        nu_p (float): proton anomalous density :math:`\\nu` [fm :sup:`-3`]
        kappa_n (float):
        kappa_p (float):

    Returns:
        float: energy functional :math:`\\epsilon_{\\pi}`
    """
    return (1./4.*(v_pi(rho_n, rho_p, 'n') + kappa_n*(np.abs(rho_grad_n))**2)*nu_n**2
            +1./4.*(v_pi(rho_n, rho_p, 'p') + kappa_p*(np.abs(rho_grad_p))**2 )* nu_p**2)

def v_pi(rho_n, rho_p, q):
    # Equation 14 from Phys Rev C 104
    """
    Calculates pairing strength :math:`\\upsilon^{pi}_q` [fm :sup:`-3`] for neutrons
    or protons, for energies below 6.5 MeV.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        q (string): nucleon type choice ('p' - proton, or 'n' - neutron)

    Returns:
        float: pairing strength :math:`\\upsilon^{pi}` [fm :sup:`-3`]

    See also:
        :func:`I`
    """
    #rho = rho_n + rho_p
    if(q=='n'):
        M = effMn(rho_n, rho_p)
    elif(q=='p'):
        M = effMp(rho_n, rho_p)
    else:
        sys.exit('# ERROR: Nucleon q must be either n or p')
    return -8.*np.pi**2/I(rho_n, rho_p, q) * (HBARC**2/2./M)**(3./2.)

def I(rho_n, rho_p, q):
    # Equation 15 from Phys Rev C 104
    """
    Approximates the solution to the integral present in the original
    :math:`\\nu^{\\pi}` formula for single-particle energies below 6.5 MeV.

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        q (string): nucleon type choice ('p' - proton, or 'n' - neutron)

    Returns:
        float: analytic integral solution

    See also:
        :func:`Lambda`
    """
    if(q=='n'):
        delta = neutron_ref_pairing_field(rho_n, rho_p)
    elif(q=='p'):
        delta = proton_ref_pairing_field(rho_n, rho_p)
    else:
        sys.exit('# ERROR: Nucleon q must be either n or p')

    mu = mu_q(rho_n, rho_p, q)
    I = mu #getting the correct size
    x = np.where(delta==NUMZERO) # intended to catch delta values with invalid kF
    I[x] = NUMZERO    
    I = np.ma.masked_where(delta == 0, I)   # masking I with incorrect delta value
                                            # (if delta=0, I = inf)
                                            # if delta(rho=0) is set to NUMZERO, I(rho=0) is too low
    y = np.where(delta>NUMZERO)
    I[y] = np.sqrt(mu[y])*(2*np.log(2*mu[y]/np.abs(delta[y]))+Lambda(6.5/mu[y]))
    return I

def Lambda(x):
    # Equation 16 from Phys Rev C 104
    """
    Function :math:`\\Lambda` used in :func:`I`.

    Args:
        x (float): variable

    Returns:
        float: `\\Lambda`

    See also:
        :func:`I`
    """
    Lambda = x
    i = np.where(x<=0 ) # intended to catch delta values with invalid kF
    Lambda[i] = NUMZERO 
    j = np.where(x>0)
    Lambda[j] = np.log(16. * x) + 2. * np.sqrt(1. + x) - 2. * np.log(1. + np.sqrt(1. + x)) - 4.
    return Lambda

def mu_q(rho_n, rho_p, q):
# Eq. taken from S. Goriely, N. Chamel, and J. M. Pearson, Phys. Rev. Lett. 102, 152503 (2009)
    """
    Calculates the chemical potential :math:`\\mu` defined with the wavevector
    :math:`k_F`.

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
