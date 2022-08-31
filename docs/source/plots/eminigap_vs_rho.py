# -*- coding: utf-8 -*-
import libnest
import numpy as np
import matplotlib.pyplot as plt
import libnest.bsk

import sys
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    sys.exit("Specify the path for the image!")

rho_n = np.linspace(0., 0.2, 10000)

e_mg = libnest.bsk.E_minigap_n(rho_n, 0.)

plt.figure()
plt.title("Energy of minigap", fontsize=15)
plt.xlabel(r"$\rho \: [{fm}^{-3}]$", fontsize=10)
plt.ylabel(r"$E_{mg} \: [MeV]$", fontsize=10)
plt.plot(rho_n, e_mg, linewidth=2.0)
plt.legend()
plt.savefig(filename)
