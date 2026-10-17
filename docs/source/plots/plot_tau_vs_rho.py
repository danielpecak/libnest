#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Kinetic density of a uniform Fermi gas, :math:`\tau = \frac{3}{5} k_F^2 \rho`, as a
function of the density."""
# Created on Tue Aug 30 14:20:50 2022 by aleks
import numpy as np
import matplotlib.pyplot as plt
from _common import output_path, savefig
import libnest.definitions


filename = output_path()

rho = np.linspace(0., 0.2, 100)

tau = libnest.definitions.rho2tau(rho)

plt.figure()
plt.title("Kinetic density", fontsize=15)
plt.xlabel(r"$\rho \: {[fm^{-3}]}$", fontsize=10)
plt.ylabel(r"$\tau \: [{fm}^{-5}]$", fontsize=10)
plt.plot(rho, tau, linewidth=2.0)
plt.xlim([0,.2])
plt.ylim([0,.4])
# plt.legend()
savefig(filename)
