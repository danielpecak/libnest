# -*- coding: utf-8 -*-
import libnest
import numpy as np
import matplotlib.pyplot as plt
from libnest.bsk import E_minigap_n

import sys
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    sys.exit("Specify the path for the image!")

rho_n = np.linspace(0., 0.08, 1000)
e_mg  = E_minigap_n(rho_n)

plt.figure()
plt.title("Energy of minigap", fontsize=15)
plt.xlabel(r"$\rho \: [{fm}^{-3}]$", fontsize=10)
plt.ylabel(r"$E_{mg} \: [MeV]$", fontsize=10)
plt.plot(rho_n, e_mg, linewidth=2.0)
plt.xlim([0,.08])
plt.ylim([0,.5])
plt.savefig(filename)
