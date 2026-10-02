#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from _common import output_path, savefig
import libnest.definitions
import libnest.bsk


filename = output_path()

kf = np.linspace(0., 1.4, 1000)
rho = libnest.definitions.kf2rho(kf)

delta_n   = libnest.bsk.neutron_pairing_field(rho)
delta_sym = libnest.bsk.symmetric_pairing_field(rho, rho)

plt.figure()
plt.title(r"Pairing field $\Delta$", fontsize=15)
plt.xlabel(r"$k_{F} \: [{fm}^{-1}]$", fontsize=10)
plt.ylabel(r"$\Delta \: [MeV]$", fontsize=10)
plt.plot(kf, delta_n, linewidth=2.0, label="NeuM")
plt.plot(kf, delta_sym, linewidth=2.0, label="SM")
plt.xlim([0,1.4])
plt.legend()
savefig(filename)
