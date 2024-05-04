#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
sys.path.insert(0, os.path.abspath('../../'))
import libnest
import numpy as np
import matplotlib.pyplot as plt
from libnest.bsk import *

import sys
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    sys.exit("Specify the path for the image!")

rho = np.linspace(0., 0.1, 1000)

C0rho = C0_rho(rho/2, rho/2)
C1rho = C1_rho(rho/2, rho/2)
C0tau = C0_tau(rho/2, rho/2)
C1tau = C1_tau(rho/2, rho/2)


plt.figure()
plt.title(r"$C^0_\rho, C^1_\rho, C^0_\tau, C^1_\tau$", fontsize=15)
plt.xlabel(r"$\rho \: [{fm}^{-3}]$", fontsize=10)
plt.ylabel(r"$C_0 \: $", fontsize=10)
plt.plot(rho, C0rho, linewidth=2.0, label=r'$C^0_\rho$')
plt.plot(rho, C1rho, linewidth=2.0, label=r'$C^1_\rho$')
plt.plot(rho, C0tau, linewidth=2.0, label=r'$C^0_\tau$')
plt.plot(rho, C1tau, linewidth=2.0, label=r'$C^1_\tau$')
plt.xlim([0,0.1])
plt.legend()
plt.savefig(filename)
