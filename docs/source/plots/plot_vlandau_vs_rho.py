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

rho = np.linspace(0., 0.1, 1000)

rho_n = rho #only for pure neutron matter
rho_p = 0

delta_n = libnest.bsk.neutron_ref_pairing_field(rho_n, rho_p)

v_landau = libnest.definitions.vLandau(delta_n, libnest.definitions.rho2kf(rho))

plt.figure()
plt.title("Landau velocity", fontsize=15)
plt.xlabel(r"$\rho \: [{fm}^{-3}]$", fontsize=10)
plt.ylabel(r"$v_{L} \: [\% \: c]$", fontsize=10)
plt.plot(rho, v_landau, linewidth=2.0)
plt.xlim([0,0.1])
plt.ylim([0,2.7])
plt.savefig(filename)
