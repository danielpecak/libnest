# -*- coding: utf-8 -*-
import libnest
import numpy as np
import matplotlib.pyplot as plt
from libnest.definitions import kf2rho
from libnest.bsk import E_minigap_rho_n

import sys
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    sys.exit("Specify the path for the image!")

kf   = np.linspace(0.001, 1.5, 1000)
e_mg = E_minigap_rho_n(kf2rho(kf))

plt.figure()
plt.title("Energy of minigap", fontsize=15)
plt.xlabel(r"$k_{F} \: [{fm}^{-1}]$", fontsize=10)
plt.ylabel(r"$E_{mg} \: [MeV]$", fontsize=10)
plt.plot(kf, e_mg, linewidth=2.0)
plt.xlim([0,1.5])
plt.ylim([0,.5])
plt.savefig(filename)
