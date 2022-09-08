# -*- coding: utf-8 -*-
import libnest
import numpy as np
import matplotlib.pyplot as plt
# import libnest.definitions

import sys
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    sys.exit("Specify the path for the image!")

sys.exit()
# kf = np.linspace(0, 1., 100)
#
# rho = libnest.definitions.kf2rho(kf)
#
# plt.figure()
# plt.title(r"Matter density against wavevector $k_F$", fontsize=15)
# plt.xlabel(r"$k_{F} \: [{fm}^{-1}]$", fontsize=10)
# plt.ylabel(r"$\rho \: [{fm}^{-3}]$", fontsize=10)
# plt.plot(kf, rho, linewidth=2.0)
# plt.legend()
# plt.savefig(filename)
