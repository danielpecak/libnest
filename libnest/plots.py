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
    or effective mass), :math:`B_{q}` [MeV fm :sup:`2`], in matter of
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
    Plots the isoscalar effective mass, :math:`M_s^*` [MeV]
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
    Plots the isovector effective mass, :math:`M_v^*`` [MeV]
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

    
def plot_epsilon(rho_n, rho_p, rho_grad, tau, j, nu, q, kappa):
    """    
    Plots the energy density :math:`\\epsilon` [MeV fm :sup:`-3`]
    in matter of density :math:`\\rho`, where :math:`\\rho` is the sum of proton
    and neutron densities, :math:`\\rho_p` and :math:`\\rho_n`.
    
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
        
    Returns:
        None
    
    See also:
        :func:`epsilon`
    """
    rho = np.linspace(0.000, 1., 100)
    rho_n = rho_n * rho
    rho_p = rho_p * rho
    rho = rho_n + rho_p
    
    grad = np.linspace(0.000, 1., 100)
    rho_grad = rho_grad * grad
    
    eps = libnest.bsk.epsilon(rho_n, rho_p, rho_grad, tau, j, nu, q, kappa)

    plt.title(r"EDF $\mathcal{E}$", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"$\mathcal{E} \:\: [MeV fm^{-3}$", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, eps, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show()
    

    
def plot_epsilon_tau(rho_n, rho_p, tau, j):
    """    
    Plots the energy density :math:`\\epsilon` [MeV fm :sup:`-3`] in neutron matter,
    related to the density-dependent effective mass.
    
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
    rho = rho_n + rho_p
    
    Mn = libnest.bsk.epsilon_tau(rho, tau, j)
    plt.title(r"E ", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"$\mathcal{E}_{\\tau} [MeV fm^{-3}$", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, Mn, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show() 
    
def plot_epsilon_delta(rho_n, rho_p, rho_grad):
    """    
    Plots the isovector effective mass, :math:`M_v^*`` [MeV]
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
    rho = rho_n + rho_p
    
    grad = np.linspace(0., 1., 100)
    rho_grad = rho_grad * grad
    
    Mn = libnest.bsk.epsilon_delta_rho(rho, rho_grad)
    plt.title(r"E", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"$\mathcal{E}_{\\Delta \\rho} [MeV fm^{-3}$", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, Mn, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show()     
    

def plot_epsilon_rho_np(rho_n, rho_p):
    """    
    Plots the isovector effective mass, :math:`M_v^*`` [MeV]
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
    rho = rho_n + rho_p
    
    Mn = libnest.bsk.g_e_rho_np(rho_n, rho_p)
    plt.title(r"E", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"$\mathcal{E}_{\\rho} [MeV fm^{-3}$", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, Mn, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show()
    
        
def plot_epsilon_tau_np(rho_n, rho_p, tau_n, tau_p, jsum2, jdiff2):
    """    
    Plots the isovector effective mass, :math:`M_v^*`` [MeV]
    in matter of density :math:`\\rho`, where :math:`\\rho` is the sum of proton
    and neutron densities, :math:`\\rho_p` and :math:`\\rho_n`.
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]
        
    Returns:
        None
    """
    rho = np.linspace(0., 1., 100)
    rho_n = rho_n * rho
    rho_p = rho_p * rho
    rho = rho_n + rho_p
    
    Mn = libnest.bsk.g_e_tau_np(rho_n, rho_p, tau_n, tau_p, jsum2, jdiff2)
    plt.title(r"E", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"$\mathcal{E}_{\\tau} [MeV fm^{-3}$", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, Mn, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show()
    
def plot_epsilon_delta_rho_np(rho_n, rho_p, rho_grad_n, rho_grad_p, rho_grad):
    """    
    Plots the isovector effective mass, :math:`M_v^*`` [MeV]
    in matter of density :math:`\\rho`, where :math:`\\rho` is the sum of proton
    and neutron densities, :math:`\\rho_p` and :math:`\\rho_n`.
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]
        
    Returns:
        None

    """
    rho = np.linspace(0., 1., 100)
    rho_n = rho_n * rho
    rho_p = rho_p * rho
    rho = rho_p + rho_n
    
    #for testing the plot
    rho_grad = rho_grad * np.linspace(0., 1., 100)
    rho_grad_n_square = (np.linspace(0., 1., 100) * rho_grad_n)**2
    rho_grad_p_square = (np.linspace(0., 1., 100) * rho_grad_p)**2
    rho_grad_square = rho_grad**2
   
    Mn = libnest.bsk.g_e_LaplaceRho_np(rho_n, rho_p, rho_grad_n_square, rho_grad_p_square, rho_grad_square)
    plt.title(r"E $\mathcal{E}$", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"$\mathcal{E}_{\\Delta \\rho} [MeV fm^{-3}]$", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, Mn, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show() 

    #to delete later
def epsilon_test(rho_n, rho_p, rho_grad_n, rho_grad_p,  tau_n, tau_p, jsum2, jdiff2):
    """    
    Plots the energy density functional from bsk_functional_full code (for testing).
    
    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
        rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
        rho_grad_n (float): neutron density gradient :math:`\\nabla \\rho` [fm :sup:`-4`]
        rho_grad_p (float): proton density gradient :math:`\\nabla \\rho` [fm :sup:`-4`]
        tau_n (float): kinetic density :math:`\\tau` [fm :sup:`-5`]
        tau_n (float): kinetic density :math:`\\tau` [fm :sup:`-5`]
        jsum2 (float): sum of momentum density/current vectors :math:`j` [fm :sup:`-3`]
        jdiff2 (float) : difference of momentum density/current vectors :math:`j` [fm :sup:`-3`]
        
    Returns:
        None

    """
    rho = np.linspace(0., 1., 100)
    rho_n = rho_n * rho
    rho_p = rho_p * rho
    rho = rho_n + rho_p
    
    rho_grad = np.linspace(0., 1., 100)
    rho_grad_n_square = (rho_grad * rho_grad_n)**2
    rho_grad_p_square = (rho_grad * rho_grad_p)**2
    rho_grad_square = rho_grad_n_square + rho_grad_p_square 
    
    x = libnest.units.HBARC**2/2/libnest.units.MN*tau_n + libnest.bsk.g_e_LaplaceRho_np(rho_n, rho_p, rho_grad_n_square, rho_grad_p_square, rho_grad_square) + libnest.bsk.g_e_tau_np(rho_n, rho_p, tau_n, tau_p, jsum2, jdiff2) + libnest.bsk.g_e_rho_np(rho_n, rho_p)
    plt.title(r"EDF$\mathcal{E}$", fontsize=15)
    plt.xlabel(r"$\rho \:\: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"$$\mathcal{E} [MeV fm^{-3}]$", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, x, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show() 

def epsilon_np(rho_n, rho_p, rho_grad_n, rho_grad_p, tau_n, tau_p, jsum2, jdiff2, nu_n, nu_p, kappa_n, kappa_p):
    """    
    Plots the energy functional :math:`\\epsilon` [MeV fm :sup:`-3`] against
    :math`\\rho`, the sum of proton and neutron of densities, :math:`\\rho_p`
    and :math:`\\rho_n` respectively.
    
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
        
    Returns:
        None
    
    See also:
        :func:`isovectorM`
    """
    rho = np.linspace(0., 1., 100)
    rho_n = rho_n * rho
    rho_p = rho_p * rho
    rho = rho_n + rho_p
    
    rho_grad = np.linspace(0., 1., 100)
    rho_grad_n = rho_grad_n * rho_grad
    rho_grad_p = rho_grad_p * rho_grad
    
    epsilon = libnest.bsk.epsilon_np(rho_n, rho_p, rho_grad_n, rho_grad_p, tau_n, tau_p, jsum2, jdiff2, nu_n, nu_p, kappa_n, kappa_p)
    plt.title(r"EDF $\mathcal{E}$", fontsize=15)
    plt.xlabel(r"$\rho \:\: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"$\mathcal{E} \:\: [MeV fm^{-3}]$", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(rho, epsilon, linewidth=2.0, label='Fit')
    #plt.legend()
    plt.show() 

    
if __name__ == '__main__':
    pass
