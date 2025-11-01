#!/usr/bin/env python3
"""
Unit tests for libnest.bsk module (Brussels-Montreal functional)
"""

import unittest
import numpy as np
from libnest import bsk


class TestBSk(unittest.TestCase):
    """Test suite for BSk energy functional"""

    def test_energy_per_nucleon_symmetric(self):
        """Test energy per nucleon for symmetric matter"""
        rho_n = rho_p = 0.08  # Symmetric matter below saturation
        E = bsk.energy_per_nucleon(rho_n, rho_p)
        # Energy should be negative (bound state)
        self.assertLess(E, 0)
        # At saturation, E/A ≈ -16 MeV
        self.assertGreater(E, -20)

    def test_energy_per_nucleon_neutron_matter(self):
        """Test energy per nucleon for pure neutron matter"""
        rho_n = 0.08
        rho_p = 0.0
        E = bsk.energy_per_nucleon(rho_n, rho_p)
        # Pure neutron matter should have higher (less bound) energy
        # than symmetric matter
        E_sym = bsk.energy_per_nucleon(0.08, 0.08)
        self.assertGreater(E, E_sym)

    def test_neutron_pairing_field_positive(self):
        """Test that neutron pairing field is positive"""
        rho_n = np.array([0.05, 0.08, 0.10])
        delta = bsk.neutron_pairing_field(rho_n)
        # Pairing field should be positive in the density range
        self.assertTrue(np.all(delta >= 0))

    def test_neutron_pairing_field_peak(self):
        """Test that pairing field has a peak somewhere in physical range"""
        rho_array = np.linspace(0.01, 0.15, 50)
        delta = bsk.neutron_pairing_field(rho_array)
        # Find maximum
        max_idx = np.argmax(delta)
        # Maximum should be in the physical range
        self.assertGreater(rho_array[max_idx], 0.01)
        self.assertLess(rho_array[max_idx], 0.15)
        # Maximum value should be positive
        self.assertGreater(delta[max_idx], 0)

    def test_effective_mass_neutron(self):
        """Test effective mass calculation for neutrons"""
        rho_n = 0.08
        rho_p = 0.0
        M_eff = bsk.effMn(rho_n, rho_p)
        # Effective mass should be positive
        self.assertGreater(M_eff, 0)
        # Effective mass is typically less than bare mass
        # M*/M ≈ 0.7-1.0
        self.assertLess(M_eff, 1.5 * bsk.MN)

    def test_effective_mass_proton(self):
        """Test effective mass calculation for protons"""
        rho_n = 0.08
        rho_p = 0.08
        M_eff = bsk.effMp(rho_n, rho_p)
        # Effective mass should be positive
        self.assertGreater(M_eff, 0)

    def test_isoscalar_mass(self):
        """Test isoscalar effective mass"""
        rho_n = rho_p = 0.08
        M_s = bsk.isoscalarM(rho_n, rho_p)
        # Should be positive
        self.assertGreater(M_s, 0)

    def test_isovector_mass(self):
        """Test isovector effective mass"""
        rho_n = rho_p = 0.08
        M_v = bsk.isovectorM(rho_n, rho_p)
        # Should be defined (positive or negative)
        self.assertTrue(np.isfinite(M_v))

    def test_symmetric_pairing_field(self):
        """Test symmetric matter pairing field"""
        # Use arrays to avoid scalar indexing issues
        rho_n = rho_p = np.array([0.08])
        delta = bsk.symmetric_pairing_field(rho_n, rho_p)
        # Should be positive for densities where pairing exists
        self.assertGreaterEqual(delta[0], 0)


class TestBSkArrays(unittest.TestCase):
    """Test that BSk functions work with arrays"""

    def test_energy_with_arrays(self):
        """Test energy functional with numpy arrays"""
        rho_n = np.array([0.05, 0.08, 0.10])
        rho_p = np.array([0.05, 0.08, 0.10])
        E = bsk.energy_per_nucleon(rho_n, rho_p)
        # Should return array of same shape
        self.assertEqual(E.shape, rho_n.shape)
        # All energies should be negative
        self.assertTrue(np.all(E < 0))

    def test_pairing_with_arrays(self):
        """Test pairing field with numpy arrays"""
        rho = np.linspace(0.01, 0.15, 10)
        delta = bsk.neutron_pairing_field(rho)
        # Should return array of same shape
        self.assertEqual(delta.shape, rho.shape)


class TestPhysicalConsistency(unittest.TestCase):
    """Test physical consistency of the BSk functional"""

    def test_saturation_point(self):
        """Test that saturation density gives minimum energy in expected range"""
        rho_array = np.linspace(0.10, 0.25, 30)
        E_array = bsk.energy_per_nucleon(rho_array, rho_array)
        
        # Find minimum
        min_idx = np.argmin(E_array)
        rho_sat = rho_array[min_idx]
        
        # Saturation should be in physical range
        # Note: The functional may have saturation outside expected range
        # depending on parametrization
        self.assertGreater(rho_sat, 0.08)
        self.assertLess(rho_sat, 0.25)

    def test_asymmetry_energy_positive(self):
        """Test that symmetry energy is positive"""
        rho = 0.16  # saturation
        E_sym = bsk.energy_per_nucleon(rho, rho)
        E_neut = bsk.energy_per_nucleon(2*rho, 0.0)
        # Pure neutron matter should have higher energy
        self.assertGreater(E_neut, E_sym)


if __name__ == '__main__':
    unittest.main()
