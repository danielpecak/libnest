# -*- coding: utf-8 -*-
"""
plots
"""

import numpy as np
import matplotlib.pyplot as plt
import libnest.bsk
import libnest.units

def plot_energy_per_nucleon(rho_n, rho_p):
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
    rho = np.linspace(0, 0.2, 100)
    rho_neum = rho #NeuM
    rho_n = 0.5 * rho
    rho_p = 0.5 * rho
    rho_sym = rho_n + rho_p #Symmetric
    
    En_neum = libnest.bsk.energy_per_nucleon(rho, 0)
    En_sym = libnest.bsk.energy_per_nucleon(rho_n, rho_p)
    energy_per_nucleon = plt.figure()
    #energy_per_nucleon.add_subplot(111)
    plt.title("Energy per nucleon", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel("E/A [MeV]", fontsize=10)
    plt.plot(rho_neum, En_neum, linewidth=2.0, label='e (NeuM)')
    plt.plot(rho_sym, En_sym, linewidth=2.0, label='e (sym)')
    plt.legend()

    plt.show()  


#plotting effective masses
def plot_effective_mass_symmetric_Mn():
    RHO = np.linspace(0., 2., 100)
    Mq = (libnest.bsk.q_effective_mass(libnest.units.MN, RHO, RHO/2)
          /(libnest.units.MP + libnest.units.MN))
    #eg if rho = rho_n
    #if p = n, rho = 2*n = 2*p, 
    EFFECTIVE_MASS_FIGURE = plt.figure()
    EFFECTIVE_MASS_FIGURE.add_subplot(111)
    plt.title("Neutron efective mass", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"M$^{*}_{n}$/M", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(RHO, Mq, linewidth=2.0, label='Fit')
    plt.legend()

    plt.show()
    
def plot_effective_mass_symmetric_Mp():
    RHO = np.linspace(0., 2., 100)
    Mq = (libnest.bsk.q_effective_mass(libnest.units.MP, RHO, RHO/2) 
          /(libnest.units.MP + libnest.units.MN))
    EFFECTIVE_MASS_FIGURE = plt.figure()
    EFFECTIVE_MASS_FIGURE.add_subplot(111)
    plt.title("Proton effective mass", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"M$^{*}_{p}$/M", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(RHO, Mq, linewidth=2.0, label='Fit')
    plt.legend()

    plt.show()

def plot_effective_mass_neutron_Mn():
    RHO = np.linspace(0., 2., 100)
    Mq = (libnest.bsk.q_effective_mass(libnest.units.MN, RHO, RHO)
              /(libnest.units.MN))
    EFFECTIVE_MASS_FIGURE = plt.figure()
    EFFECTIVE_MASS_FIGURE.add_subplot(111)
    plt.title("Neutron effective mass", fontsize=15)
    plt.xlabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
    plt.ylabel(r"M$^{*}_{p}$/M", fontsize=10)
    plt.xticks(fontsize=10)
    plt.plot(RHO, Mq, linewidth=2.0, label='Fit')
    plt.legend()

    plt.show()

    
if __name__ == '__main__':
    pass
