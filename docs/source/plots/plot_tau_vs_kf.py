#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import libnest
import numpy as np
import matplotlib.pyplot as plt
from libnest.definitions import kf2rho, rho2tau

import sys
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    sys.exit("Specify the path for the image!")


kf  = np.linspace(0, 1.5, 300)
tau = rho2tau(kf2rho(kf))

plt.figure()
plt.title("Kinetic density", fontsize=15)
plt.xlabel(r"$k_{F} \: [{fm}^{-1}]$", fontsize=10)
plt.ylabel(r"$\tau \: [{fm}^{-5}]$", fontsize=10)
plt.plot(kf, tau, linewidth=2.0)
plt.xlim([0,1.5])
# plt.ylim([0,.08])
plt.yscale('log')
plt.savefig(filename)
