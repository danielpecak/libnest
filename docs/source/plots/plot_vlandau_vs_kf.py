#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from _common import output_path, savefig
import libnest.definitions
import libnest.bsk


filename = output_path()

kf = np.linspace(1e-3, 2., 1000)

rho_n = libnest.definitions.kf2rho(kf) #only for pure neutron matter

delta_n = libnest.bsk.neutron_ref_pairing_field(rho_n, 0.)

v_landau = 100*libnest.definitions.vLandau(delta_n, kf) # in % of c

plt.figure()
plt.title("Landau velocity", fontsize=15)
plt.xlabel(r"$k_{F} \: [{fm}^{-1}]$", fontsize=10)
plt.ylabel(r"$v_{L} \: [\% \: c]$", fontsize=10)
plt.plot(kf, v_landau, linewidth=2.0)
plt.xlim([0,1.5])
plt.ylim([0,2.7])
savefig(filename)
