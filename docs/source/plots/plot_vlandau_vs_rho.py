#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from _common import output_path, savefig
import libnest.definitions
import libnest.bsk

filename = output_path()

rho = np.linspace(1e-6, 0.1, 10000)

rho_n = rho #only for pure neutron matter
rho_p = 0

delta_n = libnest.bsk.neutron_ref_pairing_field(rho_n, rho_p)

v_landau = 100*libnest.definitions.vLandau(delta_n, libnest.definitions.rho2kf(rho)) # in % of c

plt.figure()
plt.title("Landau velocity", fontsize=15)
plt.xlabel(r"$\rho \: [{fm}^{-3}]$", fontsize=10)
plt.ylabel(r"$v_{L} \: [\% \: c]$", fontsize=10)
plt.plot(rho, v_landau, linewidth=2.0)
plt.xlim([0,0.1])
plt.ylim([0,2.7])
savefig(filename)
