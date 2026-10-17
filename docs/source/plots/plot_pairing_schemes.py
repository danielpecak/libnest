#!/usr/bin/env python3
"""Reference gaps of neutrons and protons in asymmetric matter for the three interpolation
schemes (Chamel 2009 for BSk, Eq. 3 of BSkG3, Eq. 6 of BSkG4), at several proton fractions."""
import numpy as np
import matplotlib.pyplot as plt
from _common import LINESTYLES, SERIES, output_path, savefig
import libnest.bsk


filename = output_path()

SCHEMES = [("chamel2009", "Chamel 2009 (BSk)"), ("bskg3_eq3", "BSkG3, Eq. 3"),
           ("bskg4_eq6", "BSkG4, Eq. 6")]
PROTON_FRACTIONS = [0.0, 0.1, 0.3, 0.5]
rho = np.linspace(1e-4, 0.2, 800)

fig, axes = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(10, 5.5))
for col, (scheme, title) in enumerate(SCHEMES):
    functional = libnest.bsk.BSK31.with_pairing(scheme=scheme)
    axes[0, col].set_title(title, fontsize=11)
    for y_p, color, style in zip(PROTON_FRACTIONS, SERIES, LINESTYLES):
        rho_n, rho_p = (1. - y_p)*rho, y_p*rho
        axes[0, col].plot(rho, functional.neutron_ref_pairing_field(rho_n, rho_p),
                          color=color, linestyle=style, linewidth=2.0,
                          label=r"$Y_p = {}$".format(y_p))
        axes[1, col].plot(rho, functional.proton_ref_pairing_field(rho_n, rho_p),
                          color=color, linestyle=style, linewidth=2.0)
axes[0, 0].set_ylabel(r"$\Delta_n$ [MeV]")
axes[1, 0].set_ylabel(r"$\Delta_p$ [MeV]")
for ax in axes[1]:
    ax.set_xlabel(r"$\rho = \rho_n + \rho_p$ [fm$^{-3}$]")
for ax in axes.flat:
    ax.set_xlim(0., 0.2)
    ax.set_ylim(bottom=0.)
    ax.grid(alpha=0.25)
axes[0, 0].legend(frameon=False, fontsize=9)
fig.tight_layout()
savefig(filename)
