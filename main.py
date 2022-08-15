#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Info: who, where, when
# Author: Daniel Pęcak <daniel.pecak@pw.edu.pl>
# Warsaw Technical University, Université Libre de Bruxelles
# On leave: Institute of Physics, Polish Academy of Sciences, Warsaw
# April 2022, Brussels
# =========== Description
# TODO
# =========== Usage example
# $ ./test_libnest.py

import libnest
import libnest.plots
import libnest.real_data_plots
import numpy as np
import matplotlib.pyplot as plt

import sys
TXT_PATH = "C:\\Users\\aleks\\OneDrive\\Dokumenty\\libnest\\txt\\"

#print(fm3togcm3(1.0))

print("Siema")

# libnest.plots.plot_energy_per_nucleon_both()

# libnest.plots.plot_energy_per_nucleon(0.5, 0.5)

# libnest.plots.plot_pairing_field_n(1., 0)

# libnest.plots.plot_pairing_field_p(0.5, 0.5)

# libnest.plots.plot_effective_mass_n(0.5, 0.5)

# libnest.plots.plot_effective_mass_p(0.5, 0.5)

# libnest.plots.plot_B_q(0.5, 0.5, 'n')

# libnest.plots.plot_U_q(0.5, 0.5, 'n')

# libnest.plots.plot_isoscalarM(1., 0.)

# libnest.plots.plot_isovectorM(1., 0.)

# libnest.plots.plot_epsilon_test(1., 0., 0.) #rho_n, tau, nu

# libnest.plots.plot_epsilon(1., 0., 0., 0., 0., 0., 'n', 0.)
#for now the graphs start from 0.001
#to do: write exceptions for 0s in the functions (to substitue with NUMZERO or delete?)

#epsilon_test and epsilon functions give the same plot for tau = 0 and nu = 0
#they are different for tau = 1 and nu = 1
#test and original start to diverge in I() and I_test() due to the first term
#it is nullified when nu = 0 in the pairing energy density functional equation

# r = np.arange(0.0001,0.2,0.0001)
# # rho_n = r
# plt.plot(r, libnest.bsk.I_test(r),label='I test')
# plt.plot(r, libnest.bsk.I(r, 0., 'n'), label='I')
# # x=np.sqrt(libnest.bsk.mu_q(rho_n, 0., 'n'))*(2*np.log(2*libnest.bsk.mu_q(rho_n, 0., 'n')/np.abs(libnest.bsk.neutron_ref_pairing_field(rho_n, 0.))))
# # plt.plot(r, x, label='first term')
# plt.legend()

#REAL PLOTS
# libnest.real_data_plots.plot_density("N2600i_density.txt")

# libnest.real_data_plots.plot_density_contour("N2600i_density.txt")

# libnest.real_data_plots.plot_pairing_field("N2600i_delta.txt")

# libnest.real_data_plots.plot_current("N960i.3_current.txt")

#REAL PLOTS - CONT - SLICES
# libnest.real_data_plots.plot_density_slice("N216_T0.00.2_density.txt")

# libnest.real_data_plots.plot_current_slice("N216_T0.56.1_current.txt")
#current too low to be visible

libnest.real_data_plots.plot_pairing_field_slice("N216_T0.00.2_delta.txt")

libnest.real_data_plots.plot_pairing_field_slice("N216_T0.02.2_delta.txt")

libnest.real_data_plots.plot_pairing_field_slice("N216_T0.44.1_delta.txt")

# libnest.real_data_plots.plot_B_q_slice("N216_T0.00.2_density.txt")

# libnest.real_data_plots.plot_U_q_slice("N216_T0.00.2_density.txt")

# libnest.real_data_plots.plot_A_slice("N216_T0.32.1_A.txt")

# Temperature - delta
libnest.real_data_plots.plot_temperature_delta('216')

# libnest.real_data_plots.plot_temperature_delta('24000')