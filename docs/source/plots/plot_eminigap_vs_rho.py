#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
sys.path.insert(0, os.path.abspath('../../'))
import libnest
import numpy as np
import matplotlib.pyplot as plt
from libnest.definitions import E_minigap_rho_n

import sys
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    sys.exit("Specify the path for the image!")

rho_n = np.linspace(0.001, 0.08, 1000)
e_mg  = E_minigap_rho_n(rho_n)

plt.figure()
plt.title("Energy of minigap", fontsize=15)
plt.xlabel(r"$\rho \: [{fm}^{-3}]$", fontsize=10)
plt.ylabel(r"$E_{mg} \: [MeV]$", fontsize=10)
plt.plot(rho_n, e_mg, linewidth=2.0)
plt.xlim([0,.08])
plt.ylim([0,.5])
plt.savefig(filename)
