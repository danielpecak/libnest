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

kf = np.linspace(0., 1.7, 1000)

e_mg = libnest.bsk.E_minigap_n(libnest.definitions.kf2rho(kf))

plt.figure()
plt.title("Energy of minigap", fontsize=15)
plt.xlabel(r"$k_{F} \: [{fm}^{-1}]$", fontsize=10)
plt.ylabel(r"$E_{mg} \: [MeV]$", fontsize=10)
plt.plot(kf, e_mg, linewidth=2.0)
plt.legend()
plt.savefig(filename)
