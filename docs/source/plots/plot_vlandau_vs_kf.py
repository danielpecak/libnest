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

kf = np.linspace(0., 2., 1000)

rho_n = libnest.definitions.kf2rho(kf) #only for pure neutron matter

delta_n = libnest.bsk.neutron_ref_pairing_field(rho_n, 0.)

v_landau = libnest.definitions.vLandau(delta_n, kf)

plt.figure()
plt.title("Landau velocity", fontsize=15)
plt.xlabel(r"$k_{F} \: [{fm}^{-1}]$", fontsize=10)
plt.ylabel(r"$v_{L} \: [\% \: c]$", fontsize=10)
plt.plot(kf, v_landau, linewidth=2.0)
plt.legend()
plt.savefig(filename)
