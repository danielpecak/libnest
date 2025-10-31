#!/usr/bin/env python3
"""
LibNest - Library for Neutron Star Physics
===========================================

This is a simple example demonstrating basic usage of the libnest library.
For more comprehensive examples, see the examples/ directory.

Author: Daniel Pęcak <daniel.pecak@pw.edu.pl>
Warsaw Technical University, Université Libre de Bruxelles
Institute of Physics, Polish Academy of Sciences, Warsaw
"""

import libnest
import libnest.bsk
import libnest.definitions
import libnest.units
import numpy as np


def main():
    """Main function demonstrating basic libnest functionality."""
    
    print("=" * 60)
    print("LibNest - Library for Neutron Star Physics")
    print("=" * 60)
    print()
    
    # Example 1: Unit conversions
    print("Example 1: Unit Conversions")
    print("-" * 40)
    rho_fm3 = 0.16  # saturation density in fm^-3
    rho_gcm3 = libnest.units.fm3togcm3(rho_fm3)
    print(f"Saturation density: {rho_fm3:.3f} fm⁻³ = {rho_gcm3:.2e} g/cm³")
    
    temp_mev = 1.0  # temperature in MeV
    temp_k = libnest.units.MeVtoK(temp_mev)
    print(f"Temperature: {temp_mev:.1f} MeV = {temp_k:.2e} K")
    print()
    
    # Example 2: Fermi wavevector calculations
    print("Example 2: Fermi Wavevector")
    print("-" * 40)
    rho_n = 0.08  # neutron density in fm^-3
    kF = libnest.definitions.rho2kf(rho_n)
    eF = libnest.definitions.eF_n(kF)
    print(f"Neutron density: {rho_n:.3f} fm⁻³")
    print(f"Fermi wavevector: {kF:.3f} fm⁻¹")
    print(f"Fermi energy: {eF:.3f} MeV")
    print()
    
    # Example 3: BSk energy functional
    print("Example 3: Brussels-Montreal Energy Functional (BSk)")
    print("-" * 40)
    rho_n = 0.08  # neutron density
    rho_p = 0.08  # proton density (symmetric matter)
    E_per_A = libnest.bsk.energy_per_nucleon(rho_n, rho_p)
    print(f"Symmetric nuclear matter (ρₙ = ρₚ = {rho_n:.2f} fm⁻³)")
    print(f"Energy per nucleon: {E_per_A:.3f} MeV")
    print()
    
    # Example 4: Pairing field
    print("Example 4: Neutron Pairing Field")
    print("-" * 40)
    rho_n = 0.05
    delta_n = libnest.bsk.neutron_pairing_field(rho_n)
    print(f"Neutron density: {rho_n:.3f} fm⁻³")
    print(f"Pairing gap: {delta_n:.3f} MeV")
    print()
    
    print("=" * 60)
    print("For more examples, see the examples/ directory")
    print("For legacy test cases, see examples/legacy_tests.py")
    print("=" * 60)


if __name__ == '__main__':
    main()
