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
import libnest.delta_and_temperature
import libnest.wdata
import numpy as np
import matplotlib.pyplot as plt

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

# libnest.plots.epsilon_np(1., 0., 1., 0., 1., 0., 0., 0., 0.0, 0., 0., 0.) #libnest.bsk.KAPPAN, libnest.bsk.KAPPAP)
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



# r = np.arange(0.0001, 4., 0.01)
# rho_n = r
# plt.plot(r, libnest.bsk.v_pi_test_P(r, 1),label='test')
# plt.plot(r, libnest.bsk.v_pi(r, 1, 'p'),label='og')
# plt.plot(r, libnest.bsk.epsilon_pi_np(r, 0, 0, 0, 0.01, 0, libnest.bsk.KAPPAN, libnest.bsk.KAPPAP), label='og')
# plt.plot(r, libnest.bsk.epsilon_pi_test(r, 0, 0, 0, 0, 0), label='test')
# plt.plot(r, libnest.bsk.neutron_ref_pairing_field(r, 0))
# print(libnest.bsk.mu_q(r, 0, 'n'))
# plt.plot(r, libnest.bsk.neutron_ref_pairing_field(r, 0.))
# plt.plot(r, libnest.bsk.rho2kf(r))
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

#Velocities
# libnest.real_data_plots.plot_vsf("N216_T0.00.2_density.txt")
# libnest.real_data_plots.plot_vsf_nv("N216_T0.000.2_density.txt", "N216_T0.00.2_A.txt", "N216_T0.00.2_current.txt" )
# libnest.real_data_plots.plot_landau_critical_velocity("N216_T0.00.2_density.txt", "N216_T0.00.2_delta.txt" )
# libnest.real_data_plots.plot_landau_velocity("N216_T0.00.2_density.txt", "N216_T0.00.2_delta.txt")
# libnest.real_data_plots.plot_landau_velocity_temp('8000')
# #issue with kf being 0 when rho is 0


#available separately: vsf, vsf_nv+v_nv, v_landau, v_landau+v_critical
#can take fv_landau rom -60 to 60
#nv *c
#
"""
#minigap theory vs real data 
#velocities for diff temps, for diff densities
#andreiev? states fir diff temps + minigap
#temp vs delta also from low r - bigger diff, changed t crit to come from uniform matter
#pressure and sound speed
#xi cisnienie i dzwiek (speed of sound later compaer to other velocities)
#intergals numerically
#eq a3 - grad/laplace = 0, uniform, rho = rho_n = rhoq, no last line
# rozniczkowac
"""


# Temperature - delta
# libnest.delta_and_temperature.plot_temperature_delta('216')

# libnest.delta_and_temperature.plot_temperature_delta('24000')

# libnest.delta_and_temperature.plot_max_delta_temperature('8000')

# libnest.delta_and_temperature.plot_max_delta_temperature_uniform('216')


#eminigap

libnest.real_data_plots.plot_eminigap('216')
