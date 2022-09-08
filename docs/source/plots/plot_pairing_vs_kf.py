# -*- coding: utf-8 -*-
import libnest
import numpy as np
import matplotlib.pyplot as plt
import libnest.definitions
import libnest.bsk

import sys
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    sys.exit("Specify the path for the image!")

rho = np.linspace(0., 0.2, 1000)

kf = libnest.definitions.rho2kf(rho)

delta_n = libnest.bsk.neutron_pairing_field(rho)

delta_sym = libnest.bsk.symmetric_pairing_field(rho/2, rho/2)

plt.figure()
plt.title(r"Pairing field $\Delta$", fontsize=15)
plt.xlabel(r"$k_{F} : [{fm}^{-1}]$", fontsize=10)
plt.ylabel(r"$\Delta [MeV]$", fontsize=10)
plt.plot(rho, delta_n, linewidth=2.0, label="NeuM")
plt.plot(rho, delta_sym, linewidth=2.0, label="symm.")
plt.legend()
plt.savefig(filename)
