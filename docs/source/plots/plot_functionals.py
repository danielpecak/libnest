#!/usr/bin/env python3
"""Energy per nucleon in symmetric and pure neutron matter, and the neutron effective mass in
neutron matter, for every functional registered in libnest.edf."""
import numpy as np
import matplotlib.pyplot as plt
from _common import LINESTYLES, SERIES, output_path, savefig
from libnest.edf import available_functionals, get_functional


filename = output_path()

names = available_functionals()
if len(names) > len(SERIES):
    raise SystemExit("more functionals than colors: split this figure (never cycle colors)")
rho = np.linspace(1e-3, 0.3, 600)

fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
for name, color, style in zip(names, SERIES, LINESTYLES):
    f = get_functional(name)
    kwargs = dict(color=color, linestyle=style, linewidth=2.0, label=name)
    axes[0].plot(rho, f.energy_per_nucleon(rho/2, rho/2), **kwargs)
    axes[1].plot(rho, f.energy_per_nucleon(rho, 0.*rho), **kwargs)
    axes[2].plot(rho, f.effMn(rho, 0.*rho), **kwargs)
axes[0].set_title("Symmetric matter", fontsize=11)
axes[0].set_ylabel(r"$E/A$ [MeV]")
axes[1].set_title("Pure neutron matter", fontsize=11)
axes[1].set_ylabel(r"$E/A$ [MeV]")
axes[2].set_title("Pure neutron matter", fontsize=11)
axes[2].set_ylabel(r"$M^*_n/M$")
for ax in axes:
    ax.set_xlabel(r"$\rho$ [fm$^{-3}$]")
    ax.set_xlim(0., 0.3)
    ax.grid(alpha=0.25)
if len(names) > 1:
    axes[0].legend(frameon=False, fontsize=9)
else:
    fig.suptitle(names[0], fontsize=12)
fig.tight_layout()
savefig(filename)
