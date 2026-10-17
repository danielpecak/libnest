import libnest
import libnest.bsk
import libnest.definitions

# Example 1: Unit conversion
rho_fm3 = 0.16  # nuclear saturation density
rho_gcm3 = libnest.units.fm3togcm3(rho_fm3)
print(f"Density: {rho_fm3} fm⁻³ = {rho_gcm3:.2e} g/cm³")

# Example 2: Calculate Fermi energy
rho_n = 0.08  # neutron density in fm⁻³
kF = libnest.definitions.rho2kf(rho_n)
eF = libnest.definitions.eF_n(kF)
print(f"Fermi energy: {eF:.2f} MeV")

# Example 3: BSk energy functional
rho_n = rho_p = 0.08  # symmetric nuclear matter
E_per_A = libnest.bsk.energy_per_nucleon(rho_n, rho_p)
print(f"Energy per nucleon: {E_per_A:.2f} MeV")

# Example 4: Neutron pairing gap
delta_n = libnest.bsk.neutron_pairing_field(rho_n)
print(f"Pairing gap: {delta_n:.2f} MeV")
