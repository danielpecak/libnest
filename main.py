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
import libnest.wdata
import numpy as np
import matplotlib.pyplot as plt

#print(fm3togcm3(1.0))

print("Siema")


#BSK
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



#EPSILON

# libnest.plots.plot_epsilon(1., 0., 1., 1., 0., 0., 'n', 0.) #rho_n, rho_p, rho_grad, tau, j, nu, q, kappa

# libnest.plots.plot_e_tau(1., 0., 1., 0.)

# libnest.plots.plot_e_delta(1., 0., 1.)
#changed t5 term

# libnest.plots.epsilon_rho_np(1., 0.)

# libnest.plots.epsilon_tau_np(1., 0., 1., 0., 0., 0.)

# libnest.plots.epsilon_delta-rho_np(1., 0., 1., 0., 1.)
# g_e_laplace agrees with the previous delta eq without the beta term.

# libnest.plots.epsilon_test(1., 0., 1., 0., 1., 0., 0., 0.) #rho_n, rho_p, grad_rho_n, grad_rho_p,  tau_n, tau_p, jsum2, jdiff2

libnest.plots.epsilon_np(1., 0., 1., 0., 1., 0., 0., 0., 0.01, 0., 0., 0.) #libnest.bsk.KAPPAN, libnest.bsk.KAPPAP)
# rho_n, rho_p, rho_grad_n, rho_grad_p, tau_n, tau_p, jsum2, jdiff2, nu_n, nu_p, kappa_n, kappa_p

#epsilon graphs agree when rho_p = 0 and rho_grad = 0

#there is a difference between epsilon_delta_rho and g_e_laplace (from bsk_functional_full) functions
#(g_e_laplace uses an extra BETA T4 term as well as semms to have the T2 term incomplete (?) )

#to include rho_p the code requires different equations (see NesT.pdf)
#TO DO: mark current equations as neutron matter only and write new ones for rho_p
#(ie epsilon_rho, epsilon_tau, epsilon_delta_rho, epsilon_pi)


#epsilon_test and epsilon functions agree
#although there is a slight diffference in I() and I_test() due to
#neutron _ref_pairing_field and neutron_pairing_field calculations differences?



# r = np.arange(0.000, 1., 0.0001)
# rho_n = r
# # plt.plot(r, libnest.bsk.v_pi_test_P(r, 1),label='test')
# # plt.plot(r, libnest.bsk.v_pi(r, 1, 'p'),label='og')
# plt.plot(r, libnest.bsk.epsilon_pi_np(r, 0, 0, 0, 0.01, 0, libnest.bsk.KAPPAN, libnest.bsk.KAPPAP), label='og')
# plt.plot(r, libnest.bsk.epsilon_pi_test(r, 0, 0, 0, 0, 0), label='test')


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

# libnest.real_data_plots.plot_pairing_field_slice("N216_T0.00.2_delta.txt")

# libnest.real_data_plots.plot_pairing_field_slice("N216_T0.02.2_delta.txt")

# libnest.real_data_plots.plot_pairing_field_slice("N216_T0.44.1_delta.txt")

# libnest.real_data_plots.plot_B_q_slice("N216_T0.00.2_density.txt")

# libnest.real_data_plots.plot_U_q_slice("N216_T0.00.2_density.txt")

# libnest.real_data_plots.plot_A_slice("N216_T0.32.1_A.txt")


# Temperature - delta
libnest.real_data_plots.plot_temperature_delta('216')

# libnest.real_data_plots.plot_temperature_delta('24000')


