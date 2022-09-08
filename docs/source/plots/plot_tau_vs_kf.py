# -*- coding: utf-8 -*-
import libnest
import numpy as np
import matplotlib.pyplot as plt
import libnest.definitions

import sys
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    sys.exit("Specify the path for the image!")


kf = np.linspace(0, 1., 100)

tau = libnest.definitions.rho2tau(libnest.definitions.kf2rho(kf))

plt.figure()
plt.title("Kinetic density", fontsize=15)
plt.xlabel(r"$k_{F} \: [{fm}^{-1}]$", fontsize=10)
plt.ylabel(r"$\tau \: [{fm}^{-5}]$", fontsize=10)
plt.plot(kf, tau, linewidth=2.0)
plt.legend()
plt.savefig(filename)
