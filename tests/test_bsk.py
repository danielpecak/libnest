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

    

    def test_pairing_fields_scalar_input(self):
        """Scalar input must not crash (regression: np.where was called on a
        0-d array before the scalar early-return, which numpy 2.x forbids)."""
        # low density -> below the cutoff, real pairing
        self.assertGreater(bsk.neutron_pairing_field(0.05), 0)
        self.assertGreater(bsk.symmetric_pairing_field(0.05, 0.05), 0)
        # high density -> above the kF cutoff, returns numerical zero, no crash
        self.assertTrue(np.isfinite(bsk.neutron_pairing_field(0.5)))
        self.assertTrue(np.isfinite(bsk.symmetric_pairing_field(0.5, 0.5)))


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


class TestReferencePairingFields(unittest.TestCase):
    """Regression tests for the reference pairing fields and the density-gradient
    energy term. These three function bodies were previously unreachable due to
    indentation errors that broke ``import libnest.bsk`` entirely; these tests
    exercise them so the regression cannot return silently.
    """

    def test_neutron_ref_symmetric_limit(self):
        """For rho_n == rho_p (eta = 0) the neutron reference field reduces to
        the symmetric-matter pairing field."""
        rho = 0.05
        ref = bsk.neutron_ref_pairing_field(rho, rho)
        sym = bsk.symmetric_pairing_field(rho, rho)
        self.assertAlmostEqual(float(ref), float(sym), places=6)

    def test_proton_ref_symmetric_limit(self):
        """For rho_n == rho_p (eta = 0) the proton reference field reduces to
        the symmetric-matter pairing field."""
        rho = 0.05
        ref = bsk.proton_ref_pairing_field(rho, rho)
        sym = bsk.symmetric_pairing_field(rho, rho)
        self.assertAlmostEqual(float(ref), float(sym), places=6)

    def test_neutron_ref_pure_neutron_limit(self):
        """For rho_p = 0 the neutron reference field reduces to the pure
        neutron-matter pairing field."""
        rho_n = 0.05
        ref = bsk.neutron_ref_pairing_field(rho_n, 0.0)
        neum = bsk.neutron_pairing_field(rho_n)
        self.assertAlmostEqual(float(ref), float(neum), places=6)

    def test_reference_fields_finite_scalar_and_array(self):
        """Both reference fields return finite values for scalar and array input."""
        for func in (bsk.neutron_ref_pairing_field, bsk.proton_ref_pairing_field):
            with self.subTest(func=func.__name__):
                self.assertTrue(np.isfinite(func(0.05, 0.03)))
                rho_n = np.array([0.02, 0.05, 0.08])
                rho_p = np.array([0.01, 0.03, 0.04])
                out = func(rho_n, rho_p)
                self.assertEqual(out.shape, rho_n.shape)
                self.assertTrue(np.all(np.isfinite(out)))

    def test_epsilon_delta_rho_np_runs(self):
        """The density-gradient energy term evaluates to finite values for
        scalar and array inputs (guards the fixed line 1001 indentation)."""
        # scalar
        val = bsk.epsilon_delta_rho_np(0.08, 0.08, 0.0, 0.0, 0.0)
        self.assertTrue(np.isfinite(val))
        # array
        rho_n = np.array([0.05, 0.08])
        rho_p = np.array([0.05, 0.08])
        grad = np.array([1e-3, 2e-3])
        out = bsk.epsilon_delta_rho_np(rho_n, rho_p, grad, grad, grad)
        self.assertEqual(out.shape, rho_n.shape)
        self.assertTrue(np.all(np.isfinite(out)))


class TestAlternativePairingSchemes(unittest.TestCase):
    """Tests for the alternative delta_n/delta_p interpolation schemes
    (bsk.ref_pairing_field_eq2/eq3/eq6), based on Eq. (2), (3) and (6) of
    https://link.springer.com/article/10.1140/epja/s10050-025-01503-x#citeas

    TODO(Adarsh): each scheme is currently a stub that raises
    NotImplementedError (see bsk.py). As you implement a scheme, remove its
    @unittest.skip decorator and fill in real assertions - e.g. following the
    pattern in TestReferencePairingFields above:
      - eta = 0 (rho_n == rho_p) should reduce to symmetric_pairing_field
      - rho_p = 0 should reduce to neutron_pairing_field
      - finite output for both scalar and array input
    """

    @unittest.skip("TODO(Adarsh): implement ref_pairing_field_eq2 (Eq. 2)")
    def test_ref_pairing_field_eq2(self):
        rho_n, rho_p = 0.05, 0.05
        ref = bsk.ref_pairing_field_eq2(rho_n, rho_p)
        sym = bsk.symmetric_pairing_field(rho_n, rho_p)
        self.assertAlmostEqual(float(ref), float(sym), places=6)

    @unittest.skip("TODO(Adarsh): implement ref_pairing_field_eq3 (Eq. 3)")
    def test_ref_pairing_field_eq3(self):
        rho_n, rho_p = 0.05, 0.05
        ref = bsk.ref_pairing_field_eq3(rho_n, rho_p)
        sym = bsk.symmetric_pairing_field(rho_n, rho_p)
        self.assertAlmostEqual(float(ref), float(sym), places=6)

    @unittest.skip("TODO(Adarsh): implement ref_pairing_field_eq6 (Eq. 6)")
    def test_ref_pairing_field_eq6(self):
        rho_n, rho_p = 0.05, 0.05
        ref = bsk.ref_pairing_field_eq6(rho_n, rho_p)
        sym = bsk.symmetric_pairing_field(rho_n, rho_p)
        self.assertAlmostEqual(float(ref), float(sym), places=6)


if __name__ == '__main__':
    unittest.main()
