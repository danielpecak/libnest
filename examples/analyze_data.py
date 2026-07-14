#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Daniel Pęcak
"""
Example: analyze your own simulation data with libnest
======================================================

A minimal, end-to-end workflow of the kind a student would run on real data:

  1. load a 2D neutron-density map from text files (libnest.myio),
  2. compute a physical quantity at every grid point with the BSk31 functional
     (here the neutron 1S0 pairing gap Delta(rho_n) and Fermi momentum kF),
  3. summarize the result and save a figure.

The sample dataset is *generated on the fly* (see generate_sample_dataset)
rather than stored in the repository, so nothing needs to be kept in git.
Run it as-is:

    python examples/analyze_data.py

To use your own data, drop files in the same format (an ``<prefix>_info.txt``
grid file and an ``<prefix>_density.txt`` "x y rho" table) somewhere and point
the LIBNEST_DATA environment variable at that directory.
"""

import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless-safe: render to a file without a display
import matplotlib.pyplot as plt

from libnest import myio, bsk, definitions

PREFIX = "sample"


def generate_sample_dataset(data_dir, prefix=PREFIX):
    """Write a small 2D neutron-density map (a Gaussian blob on an 8x8 grid)
    into ``data_dir`` as ``<prefix>_info.txt`` and ``<prefix>_density.txt``.

    This stands in for a real WSLDA/WBSK simulation slice; it is cheap to make,
    so we generate it at runtime instead of committing data files.
    """
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    NX = NY = 8
    NZ = 1
    DX = DY = DZ = 1.0
    rho0, sigma = 0.085, 2.5          # peak density [fm^-3], blob width [fm]
    cx = cy = (NX - 1) / 2.0

    with open(data_dir / f"{prefix}_info.txt", "w") as f:
        f.write("# Grid information for the sample dataset\n")
        for k, v in [("NX", NX), ("NY", NY), ("NZ", NZ),
                     ("DX", DX), ("DY", DY), ("DZ", DZ)]:
            f.write(f"{k} = {v}\n")

    with open(data_dir / f"{prefix}_density.txt", "w") as f:
        f.write("# x [fm]   y [fm]   rho_n [fm^-3]\n")
        for i in range(NX):
            for j in range(NY):
                x, y = i * DX, j * DY
                rho = rho0 * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * sigma ** 2))
                f.write(f"{x:6.2f} {y:6.2f} {rho:12.6e}\n")
    return data_dir


def main():
    # Data location: LIBNEST_DATA if the user set it, else a generated sample.
    if "LIBNEST_DATA" in os.environ:
        data_dir = Path(os.environ["LIBNEST_DATA"])
    else:
        data_dir = generate_sample_dataset(Path(__file__).resolve().parent / "data")
        print(f"generated sample dataset in {data_dir}")

    # 1. Load grid info and the density map.
    dims = myio.readDimTxt(str(data_dir) + os.sep, PREFIX)
    print(f"grid [NX, NY, NZ] = {dims}")

    df = myio.txt2df(str(data_dir / PREFIX), "density", ["x", "y", "rho"])
    print(f"loaded {len(df)} grid points; "
          f"rho_n in [{df['rho'].min():.3e}, {df['rho'].max():.3e}] fm^-3")

    # 2. Compute libnest quantities per grid point.
    rho_n = df["rho"].to_numpy()
    df["kF"] = definitions.rho2kf(rho_n)
    df["gap"] = bsk.neutron_pairing_field(rho_n)
    peak = df.loc[df["gap"].idxmax()]
    print(f"neutron pairing gap: peak {peak['gap']:.3f} MeV "
          f"at rho_n = {peak['rho']:.3e} fm^-3 (kF = {peak['kF']:.3f} fm^-1)")

    # 3. Plot the pairing-gap map and save it.
    nx, ny = dims[0], dims[1]
    gap_map = df["gap"].to_numpy().reshape(nx, ny)
    plt.figure()
    plt.title(r"Neutron pairing gap $\Delta$ from a density map")
    plt.xlabel("x [fm]")
    plt.ylabel("y [fm]")
    im = plt.imshow(gap_map.T, origin="lower")
    plt.colorbar(im, label=r"$\Delta$ [MeV]")
    out = Path(__file__).resolve().parent / "analyze_data_gap.png"
    plt.savefig(out, dpi=120)
    print(f"saved figure: {out}")


if __name__ == "__main__":
    main()
