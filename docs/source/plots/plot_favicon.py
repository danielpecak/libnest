#!/usr/bin/env python3
"""Favicon of the documentation: the reference gap of pure neutron matter."""
import numpy as np
import matplotlib.pyplot as plt
from _common import SERIES, output_path
import libnest.bsk
import libnest.definitions


filename = output_path()

kf = np.linspace(0., 1.45, 300)
gap = libnest.bsk.neutron_pairing_field(libnest.definitions.kf2rho(kf))

fig = plt.figure(figsize=(1., 1.))
ax = fig.add_axes([0., 0., 1., 1.])
ax.fill_between(kf, gap, color=SERIES[0])
ax.set_xlim(0., 1.45)
ax.set_ylim(0., 1.05*gap.max())
ax.axis("off")
fig.savefig(filename, dpi=64, transparent=True)
plt.close(fig)
