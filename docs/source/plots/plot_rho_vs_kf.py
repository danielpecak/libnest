#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
sys.path.insert(0, os.path.abspath('../../'))

import libnest
import numpy as np
import matplotlib.pyplot as plt
import libnest.definitions

import sys
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    sys.exit("Specify the path for the image!")

kf  = np.linspace(1e-5, 2., 1000)
rho = np.log10(libnest.definitions.kf2rho(kf))

plt.figure()
plt.title("Density", fontsize=15)
plt.xlabel(r"$k_{F} \: [{fm}^{-1}]$", fontsize=10)
plt.ylabel(r"$\log(\rho) \: [{fm}^{-3}]$", fontsize=10)
plt.plot(kf, rho, linewidth=2.0)
plt.xlim([0,2])
plt.ylim([-10,0])
# plt.legend()
plt.savefig(filename)
