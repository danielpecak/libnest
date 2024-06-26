#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 14:20:50 2022

@author: aleks
"""
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

rho = np.linspace(0., 0.2, 100)

tau = libnest.definitions.rho2tau(rho)

plt.figure()
plt.title("Kinetic density", fontsize=15)
plt.xlabel(r"$\rho \: {[fm^{-3}]}$", fontsize=10)
plt.ylabel(r"$\tau \: [{fm}^{-5}]$", fontsize=10)
plt.plot(rho, tau, linewidth=2.0)
plt.xlim([0,.2])
plt.ylim([0,.2])
# plt.legend()
plt.savefig(filename)
