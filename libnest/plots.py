# -*- coding: utf-8 -*-
"""
plots
"""

import numpy as np
import matplotlib.pyplot as plt
import libnest.bsk
import libnest.units


# ================================
#         Energy/nucleon
# ================================

def plot_energy_per_nucleon(rho_n, rho_p):
    """    
    Plots the energy per nucleon for uniform matter of density :math:`\\rho`,
    the sum of proton and neutron densities, :math:`\\rho_p` and :math:`\\rho_n`
    respectively.
    
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]
        
    Returns:
        None
    
    See also:
        :func:`energy_per_nucleon`
    """
    rho = np.linspace(0., 1., 100)
    rho_n = rho * rho_n
    rho_p = rho * rho_p
    En = libnest.bsk.energy_per_nucleon(rho_n, rho_p)
    energy_per_nucleon = plt.figure()
    energy_per_nucleon.add_subplot(111)
    plt.title("Energy per nucleon", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel("E/A [MeV]", fontsize=10)
    plt.plot(rho, En, linewidth=2.0, label='Fit')
    plt.legend()

    plt.show()
    
def plot_energy_per_nucleon_both():
    """    
    Plots the energy per nucleon for symmetric uniform matter and neutron
    uniform matter of density :math:`\\rho`, which is the sum of proton and
    neutron densities, :math:`\\rho_p` and :math:`\\rho_n` respectively.
    
    Args:
        None
        
    Returns:
        None
    
    See also:
        :func:`energy_per_nucleon`
    """
    rho = np.linspace(0, 0.2, 100)
    rho_neum = rho #NeuM
    rho_n = 0.5 * rho
    rho_p = 0.5 * rho
    rho_sym = rho_n + rho_p #Symmetric
    
    En_neum = libnest.bsk.energy_per_nucleon(rho, 0)
    En_sym = libnest.bsk.energy_per_nucleon(rho_n, rho_p)
    energy_per_nucleon = plt.figure()
    energy_per_nucleon.add_subplot(111)
    plt.title("Energy per nucleon", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel("E/A [MeV]", fontsize=10)
    plt.plot(rho_neum, En_neum, linewidth=2.0, label='e (NeuM)')
    plt.plot(rho_sym, En_sym, linewidth=2.0, label='e (sym)')
    plt.legend()

    plt.show()  



# ================================
#         Pairing fields
# ================================
def plot_pairing_field_n(rho_n, rho_p):
    """    
    Plots the pairing field :math:`\\Delta` [MeV] for neutrons in matter of
    density :math:`\\rho`, the sum of proton and neutron densities,
    :math:`\\rho_p` and :math:`\\rho_n` respectively.
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]
        
    Returns:
        None
    
    See also:
        :func:`neutron_ref_pairing_field`
    """
    rho = np.linspace(0., 0.09, 100)
    rho_n = rho * rho_n
    rho_p = rho * rho_p
    delta = libnest.bsk.neutron_ref_pairing_field(rho_n, rho_p)
    pairing_field = plt.figure()
    pairing_field.add_subplot(111)
    plt.title("Pairing Field - neutrons", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel("pairing field $\Delta \: [MeV]$", fontsize=10)
    plt.plot(rho, delta, linewidth=2.0, label='Fit')

    plt.show()
    
def plot_pairing_field_p(rho_n, rho_p):
    """    
    Plots the pairing field :math:`\\Delta` [MeV] for protons in matter of
    density :math:`\\rho`, the sum of proton and neutron densities,
    :math:`\\rho_p` and :math:`\\rho_n` respectively.
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]
        
    Returns:
        None
    
    See also:
        :func:`proton_ref_pairing_field`
    """
    rho = np.linspace(0., 0.2, 100)
    rho_n = rho * rho_n
    rho_p = rho * rho_p
    delta = libnest.bsk.proton_ref_pairing_field(rho_n, rho_p)
    pairing_field = plt.figure()
    pairing_field.add_subplot(111)
    plt.title("Pairing Field - protons", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel("pairing field $\Delta \: [MeV]$", fontsize=10)
    plt.plot(rho, delta, linewidth=2.0, label='Fit')

    plt.show()

# ================================
#        Effective masses
# ================================


def plot_effective_mass_n(rho_n, rho_p):
    """    
    Plots the effective mass of neutron, :math:`M_{n}^{*} / M` [MeV] in matter of
    density :math:`\\rho`. :math:`\\rho` is the sum of proton and neutron densities,
    :math:`\\rho_p` and :math:`\\rho_n` respectively, and :math:`M` is the sum
    of neutron and proton masses.
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]
        
    Returns:
        None
    
    See also:
        :func:`effMn`
    """
    rho = np.linspace(0., 1., 100)
    rho_n = rho_n * rho
    rho_p = rho_p * rho
    Mn = libnest.bsk.effMn(rho_n, rho_p)
    plt.title("Neutron efective mass", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"M$^{*}_{n}$/M", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, Mn, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show()
    
    
def plot_effective_mass_p(rho_n, rho_p):
    """    
    Plots the effective mass of proton, :math:`M_{p}^{*}/ M` [MeV] in matter of
    density :math:`\\rho`. :math:`\\rho` is the sum of proton and neutron densities,
    :math:`\\rho_p` and :math:`\\rho_n` respectively, and :math:`M` is the sum
    of neutron and proton masses.
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]
        
    Returns:
        None
    
    See also:
        :func:`effMp`
    """
    rho = np.linspace(0., 1., 100)
    rho_n = rho_n * rho
    rho_p = rho_p * rho
    Mn = libnest.bsk.effMp(rho_n, rho_p)
    plt.title("Proton efective mass", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"M$^{*}_{p}$/M", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, Mn, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show()
    
    
def plot_B_q(rho_n, rho_p, q):
    """    
    Plots the mean field potential (from variation over kinetic density,
    or effective mass), :math:`B_{q}` [MeV fm:sup:`2`], in matter of
    density :math:`\\rho`, where :math:`\\rho` is the sum of proton and neutron
    densities, :math:`\\rho_p` and :math:`\\rho_n`.
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]
        q (string): nucleon type choice (neutron 'n' or proton 'p')
        
    Returns:
        None
    
    See also:
        :func:`B_q`
    """
    rho = np.linspace(0., 1., 100)
    rho_n = rho_n * rho
    rho_p = rho_p * rho
    Mn = libnest.bsk.B_q(rho_n, rho_p, q)
    plt.title(r"Mean potential field B$_q$", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"B$_{q}$", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, Mn, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show()
    
def plot_U_q(rho_n, rho_p, q):
    """    
    Plots the mean field potential (from variation over density :math:`\\rho`),
    :math:`U_{q}` [MeV] in matter of density :math:`\\rho`, where :math:`\\rho`
    is the sum of proton and neutron densities, :math:`\\rho_p` and :math:`\\rho_n`.
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]
        q (string): nucleon type choice (neutron 'n' or proton 'p')
        
    Returns:
        None
    
    See also:
        :func:`U_q`
    """
    rho = np.linspace(0., 1., 100)
    rho_n = rho_n * rho
    rho_p = rho_p * rho
    Mn = libnest.bsk.U_q(rho_n, rho_p, q)
    plt.title(r"Mean potential field $U_q$", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"U$_{q}$", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, Mn, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show()
    
def plot_isoscalarM(rho_n, rho_p):
    """    
    Plots the isoscalar effective mass, :math:`M :sub:`s`:sup:`*`` [MeV]
    in matter of density :math:`\\rho`, where :math:`\\rho` is the sum of proton
    and neutron densities, :math:`\\rho_p` and :math:`\\rho_n`.
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]
        
    Returns:
        None
    
    See also:
        :func:`isoscalarM`
    """
    rho = np.linspace(0., 1., 100)
    rho_n = rho_n * rho
    rho_p = rho_p * rho
    Mn = libnest.bsk.isoscalarM(rho_n, rho_p)
    plt.title(r"Effective isoscalar mass $M^*_s$", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"M$^*_{s}/M$", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, Mn, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show()
    
def plot_isovectorM(rho_n, rho_p):
    """    
    Plots the isovector effective mass, :math:`M :sub:`s`:sup:`*`` [MeV]
    in matter of density :math:`\\rho`, where :math:`\\rho` is the sum of proton
    and neutron densities, :math:`\\rho_p` and :math:`\\rho_n`.
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]
        
    Returns:
        None
    
    See also:
        :func:`isovectorM`
    """
    rho = np.linspace(0., 1., 100)
    rho_n = rho_n * rho
    rho_p = rho_p * rho
    Mn = libnest.bsk.isovectorM(rho_n, rho_p)
    plt.title(r"Effective isovector mass $M^*_v$", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"M$^*_{v}/M$", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, Mn, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show()
    
    
# def plot_effective_mass_symmetric_Mn():
#     RHO = np.linspace(0., 2., 100)
#     Mq = (libnest.bsk.q_effective_mass(libnest.units.MN, RHO, RHO/2)
#           /(libnest.units.MP + libnest.units.MN))
#     #eg if rho = rho_n
#     #if p = n, rho = 2*n = 2*p, 
#     EFFECTIVE_MASS_FIGURE = plt.figure()
#     EFFECTIVE_MASS_FIGURE.add_subplot(111)
#     plt.title("Neutron efective mass", fontsize=15)
#     plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
#     plt.ylabel(r"M$^{*}_{n}$/M", fontsize=10)
#     plt.xticks(fontsize=10)
#     plt.plot(RHO, Mq, linewidth=2.0, label='Fit')
#     plt.legend()

#     plt.show()
    
# def plot_effective_mass_symmetric_Mp():
#     RHO = np.linspace(0., 2., 100)
#     Mq = (libnest.bsk.q_effective_mass(libnest.units.MP, RHO, RHO/2) 
#           /(libnest.units.MP + libnest.units.MN))
#     EFFECTIVE_MASS_FIGURE = plt.figure()
#     EFFECTIVE_MASS_FIGURE.add_subplot(111)
#     plt.title("Proton effective mass", fontsize=15)
#     plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
#     plt.ylabel(r"M$^{*}_{p}$/M", fontsize=10)
#     plt.xticks(fontsize=10)
#     plt.plot(RHO, Mq, linewidth=2.0, label='Fit')
#     plt.legend()

#     plt.show()

# def plot_effective_mass_neutron_Mn():
#     RHO = np.linspace(0., 2., 100)
#     Mq = (libnest.bsk.q_effective_mass(libnest.units.MN, RHO, RHO)
#               /(libnest.units.MN))
#     EFFECTIVE_MASS_FIGURE = plt.figure()
#     EFFECTIVE_MASS_FIGURE.add_subplot(111)
#     plt.title("Neutron effective mass", fontsize=15)
#     plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
#     plt.ylabel(r"M$^{*}_{p}$/M", fontsize=10)
#     plt.xticks(fontsize=10)
#     plt.plot(RHO, Mq, linewidth=2.0, label='Fit')
#     plt.legend()

#     plt.show()

    
if __name__ == '__main__':
    pass
