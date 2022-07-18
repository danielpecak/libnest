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

#print(fm3togcm3(1.0))

print("Siema")

#print(dir(libnest))

#libnest.plots.plot_energy_per_nucleon_both()

#libnest.plots.plot_energy_per_nucleon()

libnest.plots.plot_pairing_field_n(1., 0)

#libnest.plots.plot_pairing_field_p(1., 0)

#libnest.plots.plot_effective_mass_symmetric_Mn() #symmetric

#libnest.plots.plot_effective_mass_symmetric_Mp() #symmetric

#libnest.plots.plot_effective_mass_neutron_Mn() #neutron matter

#libnest.real_data_plots.density("N2600i_density.txt")