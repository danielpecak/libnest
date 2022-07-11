# -*- coding: utf-8 -*-
"""
plots
"""

import numpy as np
import matplotlib.pyplot as plt
import libnest.definitions

#plotting effective masses

def plot_effective_mass_symmetric_Mn():
    RHO = np.linspace(0., 2., 100)
    Mq = (libnest.definitions.q_effective_mass(libnest.units.MN, RHO, RHO/2)
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
    Mq = (libnest.definitions.q_effective_mass(libnest.units.MP, RHO, RHO/2) 
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
    Mq = (libnest.definitions.q_effective_mass(libnest.units.MN, RHO, RHO)
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
