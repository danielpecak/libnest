#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys
import os
sys.path.insert(0, os.path.abspath('../../'))
import libnest
import numpy as np
import matplotlib.pyplot as plt
import libnest.bsk


import sys
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    sys.exit("Specify the path for the image!")

rho = np.linspace(0, 1., 200)
rho_neum = rho #NeuM
rho_n = 0.5 * rho
rho_p = 0.5 * rho
rho_sym = rho_n + rho_p #Symmetric

En_neum = libnest.bsk.energy_per_nucleon(rho, 0)
En_sym = libnest.bsk.energy_per_nucleon(rho_n, rho_p)
energy_per_nucleon = plt.figure()
energy_per_nucleon.add_subplot(111)
plt.title("Energy per nucleon (BSk31)", fontsize=15)
plt.xlabel(r"$\rho \: {[fm^{-3}]}$", fontsize=10)
plt.ylabel("E/A [MeV]", fontsize=10)
plt.plot(rho_neum, En_neum, linewidth=2.0, label='NeuM')
plt.plot(rho_sym, En_sym, linewidth=2.0, label='SM')
plt.xlim([0,1])
plt.legend()
plt.savefig(filename)
