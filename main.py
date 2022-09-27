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
# import libnest.wdata
import numpy as np
import matplotlib.pyplot as plt

print("Siema")


#BSK
# libnest.plots.plot_energy_per_nucleon_both()
# libnest.plots.plot_energy_per_nucleon(0.5, 0.5)


# libnest.plots.plot_pairing_field_n(0.5, 0.5)
# libnest.plots.plot_pairing_field_p(0.5, 0.5)
# libnest.plots.plot_effective_mass_n(0.5, 0.5)
# libnest.plots.plot_effective_mass_p(0.5, 0.5)

# libnest.plots.plot_pairing_field_n(1., 0)
# libnest.plots.plot_pairing_field_p(1., 0.)
# libnest.plots.plot_effective_mass_n(1., 0.)


# libnest.plots.plot_B_q(1., 0., 'n')
# libnest.plots.plot_U_q(1., 0., 'n')
# libnest.plots.plot_isoscalarM(1., 0.)
# libnest.plots.plot_isovectorM(1., 0.)



#EPSILON
# libnest.plots.plot_epsilon(1., 0., 0., 0., 0., 0., 'n', 0.) #rho_n, rho_p, rho_grad, tau, j, nu, q, kappa

# libnest.plots.epsilon_rho_np(1., 0.)
# libnest.plots.epsilon_tau_np(1., 0., 1., 0., 0., 0.)
# libnest.plots.epsilon_delta_rho_np(1., 0., 1., 0., 1.)
# g_e_laplace agrees with the previous delta eq without the beta term.

# libnest.plots.epsilon_np(1., 0., 0., 0., 0., 0., 0., 0., 0.0, 0., 0., 0.)
# libnest.plots.epsilon_np(1., 0., 1., 0., 0., 0., 0., 0., 0.0, 0., 0., 0.)
# libnest.plots.epsilon_np(1., 0., 1., 0., 1., 0., 0., 0., 0.0, 0., 0., 0.)
# rho_n, rho_p, rho_grad_n, rho_grad_p, tau_n, tau_p, jsum2, jdiff2, nu_n, nu_p, kappa_n, kappa_p

#epsilon graphs agree when rho_p = 0 and rho_grad = 0



#REAL PLOTS
# libnest.real_data_plots.plot_density("N2600i_density.txt")
# libnest.real_data_plots.plot_density_contour("N2600i_density.txt")
# libnest.real_data_plots.plot_pairing_field("N2600i_delta.txt")
# libnest.real_data_plots.plot_current("N2600i.3_current.txt")


#REAL PLOTS - CONT - SLICES
# libnest.real_data_plots.plot_density_slice("N24000_T0.36.1_density.txt")
# libnest.real_data_plots.plot_density_slice("N24000_T0.00.3_density.txt")
# libnest.real_data_plots.plot_density_slice("N24000_T0.60.1_density.txt")

# libnest.real_data_plots.plot_density_slice("N216_T0.00.2_density.txt")
# libnest.real_data_plots.plot_density_slice("N216_T0.36.1_density.txt")
# libnest.real_data_plots.plot_density_slice("N216_T0.60.1_density.txt")



# libnest.real_data_plots.plot_pairing_field_slice("N24000_T0.36.1_delta.txt")
# libnest.real_data_plots.plot_pairing_field_slice("N24000_T0.00.3_delta.txt")
# libnest.real_data_plots.plot_pairing_field_slice("N24000_T0.60.1_delta.txt")

# libnest.real_data_plots.plot_pairing_field_slice("N216_T0.00.2_delta.txt")
# libnest.real_data_plots.plot_pairing_field_slice("N216_T0.36.1_delta.txt")
# libnest.real_data_plots.plot_pairing_field_slice("N216_T0.60.1_delta.txt")



# libnest.real_data_plots.plot_B_q_slice("N24000_T0.36.1_density.txt")
# libnest.real_data_plots.plot_B_q_slice("N24000_T0.00.3_density.txt")
# libnest.real_data_plots.plot_B_q_slice("N24000_T0.60.1_density.txt")

# libnest.real_data_plots.plot_U_q_slice("N24000_T0.36.1_density.txt")
# libnest.real_data_plots.plot_U_q_slice("N24000_T0.00.3_density.txt")
# libnest.real_data_plots.plot_U_q_slice("N24000_T0.60.1_density.txt")

# libnest.real_data_plots.plot_A_slice("N24000_T0.36.1_A.txt")
# libnest.real_data_plots.plot_A_slice("N24000_T0.00.3_A.txt")
# libnest.real_data_plots.plot_A_slice("N24000_T0.60.1_A.txt")


# libnest.real_data_plots.plot_B_q_slice("N216_T0.36.1_density.txt")
# libnest.real_data_plots.plot_B_q_slice("N216_T0.00.2_density.txt")
# libnest.real_data_plots.plot_B_q_slice("N216_T0.60.1_density.txt")

# libnest.real_data_plots.plot_U_q_slice("N216_T0.36.1_density.txt")
# libnest.real_data_plots.plot_U_q_slice("N216_T0.00.2_density.txt")
# libnest.real_data_plots.plot_U_q_slice("N216_T0.60.1_density.txt")

# libnest.real_data_plots.plot_A_slice("N216_T0.36.1_A.txt")
# libnest.real_data_plots.plot_A_slice("N216_T0.00.2_A.txt")
# libnest.real_data_plots.plot_A_slice("N216_T0.60.1_A.txt")



# libnest.real_data_plots.plot_current_slice("N216_T0.56.1_current.txt")
# libnest.real_data_plots.plot_current_slice("N216_T0.36.1_current.txt")
# libnest.real_data_plots.plot_current_slice("N216_T0.00.2_current.txt")

# libnest.real_data_plots.plot_current_slice("N24000_T0.60.1_current.txt")
# libnest.real_data_plots.plot_current_slice("N24000_T0.36.1_current.txt")
# libnest.real_data_plots.plot_current_slice("N24000_T0.00.3_current.txt")
#current too low to be visible



#Velocities
# libnest.real_data_plots.plot_vsf("N216_T0.00.2_density.txt")
# libnest.real_data_plots.plot_vsf_nv("N216_T0.000.2_density.txt", "N216_T0.00.2_A.txt", "N216_T0.00.2_current.txt" )
# libnest.real_data_plots.plot_landau_critical_velocity("N24000_T0.00.3_density.txt", "N24000_T0.00.3_delta.txt" )
# libnest.real_data_plots.plot_landau_critical_velocity("N24000_T0.36.1_density.txt", "N24000_T0.36.1_delta.txt" )
# libnest.real_data_plots.plot_landau_critical_velocity("N24000_T0.60.1_density.txt", "N24000_T0.60.1_delta.txt" )

# libnest.real_data_plots.plot_landau_velocity("N216_T0.00.2_density.txt", "N216_T0.00.2_delta.txt")
# libnest.real_data_plots.plot_landau_velocity_temperature('216')
# libnest.real_data_plots.plot_landau_velocity_temperature('8000')
# libnest.real_data_plots.plot_landau_velocity_temperature('24000')
# libnest.real_data_plots.plot_speed_of_sound("N24000_T0.60.1_density.txt")
# #issue with kf being 0 when rho is 0

# libnest.plots.plot_speed_of_sound_n(1.1)
# libnest.plots.plot_speed_of_sound_n(0.09)
# libnest.plots.plot_pressure_n(0.09)
# libnest.plots.plot_v_landau(0.09)
# libnest.plots.plot_v_critical(0.09)
# libnest.plots.plot_v_sf(1.)



# Temperature - delta
# libnest.delta_and_temperature.plot_temperature_delta('24000')
# libnest.delta_and_temperature.plot_max_delta_temperature('24000')
# libnest.delta_and_temperature.plot_max_delta_temperature_uniform('24000')



#eminigap
# libnest.real_data_plots.plot_e_minigap_temperature('216')
# libnest.real_data_plots.plot_e_minigap_temperature('4160')
# libnest.real_data_plots.plot_e_minigap_temperature('8000')
# libnest.real_data_plots.plot_e_minigap_temperature('13600')
# libnest.real_data_plots.plot_e_minigap_temperature('16640')
# libnest.real_data_plots.plot_e_minigap_temperature('24000')

# libnest.real_data_plots.andreev_e_minimum("N15000i_states.0000.txt")
