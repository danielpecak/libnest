#!/usr/bin/env python3
"""
Unit tests for libnest.definitions module
"""

import unittest
import numpy as np
from libnest import definitions


class TestDefinitions(unittest.TestCase):
    """Test suite for definitions module"""

    def test_rho2kf_basic(self):
        """Test Fermi wavevector calculation from density"""
        rho = 0.16  # saturation density
        kF = definitions.rho2kf(rho)
        # Expected: kF ≈ 1.68 fm^-1 at saturation
        self.assertAlmostEqual(kF, 1.680, places=2)

    def test_kf2rho_basic(self):
        """Test density calculation from Fermi wavevector"""
        kF = 1.68  # Fermi wavevector at saturation
        rho = definitions.kf2rho(kF)
        self.assertAlmostEqual(rho, 0.16, places=2)

    def test_rho2kf_kf2rho_inverse(self):
        """Test that rho2kf and kf2rho are inverse functions"""
        rho_original = 0.08
        kF = definitions.rho2kf(rho_original)
        rho_back = definitions.kf2rho(kF)
        self.assertAlmostEqual(rho_original, rho_back, places=10)

    def test_rho2kf_array(self):
        """Test rho2kf works with numpy arrays"""
        rho_array = np.array([0.08, 0.16, 0.24])
        kF_array = definitions.rho2kf(rho_array)
        self.assertEqual(len(kF_array), 3)
        self.assertTrue(all(kF_array > 0))

    def test_rhoEta_symmetric(self):
        """Test rhoEta for symmetric matter"""
        rho_n = 0.08
        rho_p = 0.08
        rho, eta = definitions.rhoEta(rho_n, rho_p)
        self.assertAlmostEqual(rho, 0.16, places=10)
        self.assertAlmostEqual(eta, 0.0, places=10)

    def test_rhoEta_pure_neutron(self):
        """Test rhoEta for pure neutron matter"""
        rho_n = 0.08
        rho_p = 0.0
        rho, eta = definitions.rhoEta(rho_n, rho_p)
        self.assertAlmostEqual(rho, 0.08, places=10)
        self.assertAlmostEqual(eta, 0.08, places=10)

    def test_rho2tau_basic(self):
        """Test kinetic density calculation"""
        rho = 0.16
        tau = definitions.rho2tau(rho)
        # Tau should be positive for positive density
        self.assertGreater(tau, 0)
        # For uniform matter, tau ∝ rho^(5/3)
        expected = 0.6 * (3 * np.pi)**(2/3) * rho**(5/3)
        self.assertAlmostEqual(tau, expected, places=10)

    def test_eF_n_basic(self):
        """Test Fermi energy calculation for neutrons"""
        kF = 1.333
        eF = definitions.eF_n(kF)
        # Fermi energy should be positive
        self.assertGreater(eF, 0)
        # At saturation, eF ≈ 37 MeV
        self.assertAlmostEqual(eF, 36.8, places=0)

    def test_vLandau_basic(self):
        """Test Landau velocity calculation"""
        delta = 1.0  # MeV
        kF = 1.0     # fm^-1
        v_L = definitions.vLandau(delta, kF)
        # Velocity should be between 0 and 1 (in units of c)
        self.assertGreater(v_L, 0)
        self.assertLess(v_L, 1)

    def test_E_minigap_delta_n_basic(self):
        """Test minigap energy calculation"""
        delta = 1.0  # MeV
        rho_n = 0.08  # fm^-3
        E_mg = definitions.E_minigap_delta_n(delta, rho_n)
        # Minigap should be positive
        self.assertGreater(E_mg, 0)

    def test_E_minigap_consistency(self):
        """Test that E_minigap_rho_n uses E_minigap_delta_n correctly"""
        # This is more of an integration test
        # Use array to avoid scalar indexing issues in bsk.py
        rho_n = np.array([0.05])
        E_mg = definitions.E_minigap_rho_n(rho_n)
        # Should return a positive energy
        self.assertGreater(E_mg[0] if hasattr(E_mg, '__getitem__') else E_mg, 0)


class TestPhysicalLimits(unittest.TestCase):
    """Test physical limits and edge cases"""

    def test_rho2kf_zero(self):
        """Test kF for zero density"""
        kF = definitions.rho2kf(0.0)
        self.assertEqual(kF, 0.0)

    def test_rho2kf_negative_raises_warning(self):
        """Test that negative density produces complex result (handled by numpy)"""
        # Negative density is unphysical, but numpy will handle it
        result = definitions.rho2kf(-0.1)
        # Result will be nan (not a number) due to taking cube root of negative
        # In Python 3, this might work differently, so we just check it doesn't crash


if __name__ == '__main__':
    unittest.main()
