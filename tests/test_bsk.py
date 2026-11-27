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
    """

    def test_neutron_ref_pairing_field_eq2(self):
        # Symmetric limit
        rho_n, rho_p = 0.05, 0.05
        ref = bsk.neutron_ref_pairing_field_eq2(rho_n, rho_p)
        sym = bsk.symmetric_pairing_field(rho_n, rho_p)
        self.assertAlmostEqual(float(ref), float(sym), places=6)
        # Pure neutron matter limit
        rho_n, rho_p = 0.05, 0.00
        ref  = bsk.neutron_ref_pairing_field_eq2(rho_n, rho_p)
        neum = bsk.neutron_pairing_field(rho_n)
        self.assertAlmostEqual(float(ref), float(neum), places=6)
        # Pure proton matter limit
        rho_n, rho_p = 0.00, 0.05
        ref   = bsk.neutron_ref_pairing_field_eq2(rho_n, rho_p)
        protm = bsk.neutron_pairing_field(rho_n)
        self.assertAlmostEqual(float(ref), float(protm), places=6)
        
    def test_proton_ref_pairing_field_eq2(self):
        # Symmetric limit
        rho_n, rho_p = 0.05, 0.05
        ref = bsk.proton_ref_pairing_field_eq2(rho_n, rho_p)
        sym = bsk.symmetric_pairing_field(rho_n, rho_p)
        self.assertAlmostEqual(float(ref), float(sym), places=6)
        # Pure neutron matter limit
        rho_n, rho_p = 0.05, 0.00
        ref   = bsk.proton_ref_pairing_field_eq2(rho_n, rho_p)
        neum = bsk.neutron_pairing_field(rho_p) # Note: this is correct, the symmetry is assumed in this model
        self.assertAlmostEqual(float(ref), float(neum), places=6)
        # Pure proton matter limit
        rho_n, rho_p = 0.00, 0.05
        ref   = bsk.proton_ref_pairing_field_eq2(rho_n, rho_p)
        protm = bsk.neutron_pairing_field(rho_p) # Note: this is correct, the symmetry is assumed in this model
        self.assertAlmostEqual(float(ref), float(protm), places=6)


    def test_neutron_ref_pairing_field_eq3(self):
        # Symmetric limit
        rho_n, rho_p = 0.05, 0.05
        ref = bsk.neutron_ref_pairing_field_eq3(rho_n, rho_p)
        sym = bsk.symmetric_pairing_field(rho_n, rho_p)
        self.assertAlmostEqual(float(ref), float(sym), places=6)
        # Pure neutron matter limit
        rho_n, rho_p = 0.05, 0.00
        ref = bsk.neutron_ref_pairing_field_eq3(rho_n, rho_p)
        neum = bsk.neutron_pairing_field(rho_n)
        self.assertAlmostEqual(float(ref), float(neum), places=6)
        # Pure proton matter limit
        rho_n, rho_p = 0.00, 0.05
        ref = bsk.neutron_ref_pairing_field_eq3(rho_n, rho_p)
        protm = bsk.neutron_pairing_field(rho_n) # Note: this is correct, the symmetry is assumed in this model
        self.assertAlmostEqual(float(ref), float(protm), places=6)
        
    def test_proton_ref_pairing_field_eq3(self):
        # Symmetric limit
        rho_n, rho_p = 0.05, 0.05
        ref = bsk.proton_ref_pairing_field_eq3(rho_n, rho_p)
        sym = bsk.symmetric_pairing_field(rho_n, rho_p)
        self.assertAlmostEqual(float(ref), float(sym), places=6)
        # Pure neutron matter limit
        rho_n, rho_p = 0.05, 0.00
        ref = bsk.proton_ref_pairing_field_eq3(rho_n, rho_p)
        neum = bsk.neutron_pairing_field(rho_p) # Note: this is correct, the symmetry is assumed in this model
        self.assertAlmostEqual(float(ref), float(neum), places=6)
        # Pure proton matter limit
        rho_n, rho_p = 0.00, 0.05
        ref = bsk.proton_ref_pairing_field_eq3(rho_n, rho_p)
        protm = bsk.neutron_pairing_field(rho_p) # Note: this is correct, the symmetry is assumed in this model
        self.assertAlmostEqual(float(ref), float(protm), places=6)

    def test_neutron_ref_pairing_field_eq6(self):
        # Symmetric limit
        rho_n, rho_p = 0.05, 0.05
        ref = bsk.neutron_ref_pairing_field_eq6(rho_n, rho_p)
        sym = bsk.symmetric_pairing_field(rho_n, rho_p)
        self.assertAlmostEqual(float(ref), float(sym), places=6)
        # Pure neutron matter limit
        rho_n, rho_p = 0.05, 0.00
        ref = bsk.neutron_ref_pairing_field_eq6(rho_n, rho_p)
        neum = bsk.neutron_pairing_field(rho_n)
        self.assertAlmostEqual(float(ref), float(neum), places=6)
        # Pure proton matter limit
        rho_n, rho_p = 0.00, 0.05
        ref = bsk.neutron_ref_pairing_field_eq6(rho_n, rho_p)
        protm = bsk.neutron_pairing_field(rho_n) # Note: this is correct, the symmetry is assumed in this model
        self.assertAlmostEqual(float(ref), float(protm), places=6)
        
    def test_proton_ref_pairing_field_eq6(self):
        # Symmetric limit
        rho_n, rho_p = 0.05, 0.05
        ref = bsk.proton_ref_pairing_field_eq6(rho_n, rho_p)
        sym = bsk.symmetric_pairing_field(rho_n, rho_p)
        self.assertAlmostEqual(float(ref), float(sym), places=6)
        # Pure neutron matter limit
        rho_n, rho_p = 0.05, 0.00
        ref = bsk.proton_ref_pairing_field_eq6(rho_n, rho_p)
        neum = bsk.neutron_pairing_field(rho_p) # Note: this is correct, the symmetry is assumed in this model
        self.assertAlmostEqual(float(ref), float(neum), places=6)
        # Pure proton matter limit
        rho_n, rho_p = 0.00, 0.05
        ref = bsk.proton_ref_pairing_field_eq6(rho_n, rho_p)
        protm = bsk.neutron_pairing_field(rho_p) # Note: this is correct, the symmetry is assumed in this model
        self.assertAlmostEqual(float(ref), float(protm), places=6)
        
    def test_reference_fields_finite_scalar_and_array(self):
        """All reference fields return finite values for scalar and array input."""
        for func in (bsk.neutron_ref_pairing_field_eq2, bsk.neutron_ref_pairing_field_eq3,bsk.neutron_ref_pairing_field_eq6,
                     bsk.proton_ref_pairing_field_eq2,bsk.proton_ref_pairing_field_eq3,bsk.proton_ref_pairing_field_eq6):
            with self.subTest(func=func.__name__):
                self.assertTrue(np.isfinite(func(0.05, 0.03)))
                rho_n = np.array([0.02, 0.05, 0.08])
                rho_p = np.array([0.01, 0.03, 0.04])
                out = func(rho_n, rho_p)
                self.assertEqual(out.shape, rho_n.shape)
                self.assertTrue(np.all(np.isfinite(out)))


class TestMeanField(unittest.TestCase):
    """U_q = U_rho_q + U_tau_q must be the derivative of epsilon_rho + epsilon_tau over rho_q
    (at fixed tau, j), and equal U_rho_q + U_tau_q of W-BSK (hpc-engine/bsk/bsk_functional.c)."""

    DENSITIES = [(0.08, 0.08), (0.04, 0.04), (0.16, 0.0), (0.10, 0.06), (0.02, 0.09)]
    # rho_n, rho_p, tau_n, tau_p, (j_n + j_p)^2, (j_n - j_p)^2
    STATES = [(0.08, 0.08, 0.10, 0.09, 0.0, 0.0), (0.16, 0.0, 0.30, 0.0, 0.0, 0.0),
              (0.10, 0.06, 0.15, 0.07, 1e-3, 4e-4), (0.02, 0.09, 0.01, 0.12, 2e-4, 1e-4),
              (0.04, 0.04, 0.05, 0.06, 2e-4, 0.0)]
    # U_rho_n + U_tau_n, U_rho_p + U_tau_p of W-BSK (BSk31, NUCLEAR_MATTER) for STATES[:4], computed
    # with bsk_functional.c through ctypes (2026-11-27). W-BSK regularizes rho with DENSEPSILON =
    # 1e-8 inside the powers, libnest with 1e-12: the two differ by ~1e-5 MeV here; with
    # DENSEPSILON = 0 on both sides they agree to 1e-15.
    WBSK = [(-59.3548745901, -58.9037822844), (-23.1525073516, -101.3032509001),
            (-48.9571252410, -68.2692649752), (-75.5839195086, -32.3518938959)]

    def check_U_rho_q_is_derivative(self, mod):
        h = 1e-6
        for rho_n, rho_p in self.DENSITIES:
            dn = (mod.epsilon_rho_np(rho_n + h, rho_p) - mod.epsilon_rho_np(rho_n - h, rho_p)) / (2 * h)
            dp = (mod.epsilon_rho_np(rho_n, rho_p + h) - mod.epsilon_rho_np(rho_n, rho_p - h)) / (2 * h)
            with self.subTest(module=mod.__name__, rho_n=rho_n, rho_p=rho_p):
                self.assertAlmostEqual(float(mod.U_rho_q(rho_n, rho_p, 'n')), float(dn), delta=1e-4)
                self.assertAlmostEqual(float(mod.U_rho_q(rho_n, rho_p, 'p')), float(dp), delta=1e-4)

    def check_U_q_is_derivative(self, mod):
        h = 1e-6
        for rho_n, rho_p, tau_n, tau_p, jsum2, jdiff2 in self.STATES:
            def eps(rn, rp):
                return (mod.epsilon_rho_np(rn, rp)
                        + mod.epsilon_tau_np(rn, rp, tau_n, tau_p, jsum2, jdiff2))
            dn = (eps(rho_n + h, rho_p) - eps(rho_n - h, rho_p)) / (2 * h)
            dp = (eps(rho_n, rho_p + h) - eps(rho_n, rho_p - h)) / (2 * h)
            args = (tau_n, tau_p, jsum2, jdiff2)
            with self.subTest(module=mod.__name__, rho_n=rho_n, rho_p=rho_p):
                self.assertAlmostEqual(float(mod.U_q(rho_n, rho_p, 'n', *args)), float(dn), delta=1e-4)
                self.assertAlmostEqual(float(mod.U_q(rho_n, rho_p, 'p', *args)), float(dp), delta=1e-4)

    def check_U_q_wbsk(self, mod):
        for (rho_n, rho_p, *args), (u_n, u_p) in zip(self.STATES, self.WBSK):
            with self.subTest(module=mod.__name__, rho_n=rho_n, rho_p=rho_p):
                self.assertAlmostEqual(float(mod.U_q(rho_n, rho_p, 'n', *args)), u_n, delta=1e-4)
                self.assertAlmostEqual(float(mod.U_q(rho_n, rho_p, 'p', *args)), u_p, delta=1e-4)

    def test_U_rho_q_bsk(self):
        self.check_U_rho_q_is_derivative(bsk)

    def test_U_rho_q_bskg(self):
        from libnest import bskg
        self.check_U_rho_q_is_derivative(bskg)

    def test_U_q_bsk(self):
        self.check_U_q_is_derivative(bsk)
        self.check_U_q_wbsk(bsk)

    def test_U_q_bskg(self):
        from libnest import bskg
        self.check_U_q_is_derivative(bskg)
        self.check_U_q_wbsk(bskg)

    def test_U_q_without_tau_is_U_rho_q(self):
        """Every term of U_tau_q is proportional to tau or j^2."""
        for rho_n, rho_p in self.DENSITIES:
            for q in 'np':
                self.assertEqual(float(bsk.U_q(rho_n, rho_p, q, 0., 0.)), float(bsk.U_rho_q(rho_n, rho_p, q)))

    def test_U_q_default_tau_is_fermi_gas(self):
        from libnest.definitions import rho2tau
        rho_n = np.array([0.0, 0.02, 0.08, 0.16])
        rho_p = np.array([0.0, 0.01, 0.08, 0.0])
        for q in 'np':
            np.testing.assert_array_equal(bsk.U_q(rho_n, rho_p, q),
                                          bsk.U_q(rho_n, rho_p, q, rho2tau(rho_n), rho2tau(rho_p)))
        self.assertTrue(np.all(np.isfinite(bsk.U_q(rho_n, rho_p, 'n'))))

    def test_U_q_symmetric_matter_bound(self):
        """Symmetric matter at saturation: U_q ~ -60 MeV (U_rho -65.8 + U_tau +6.0), n and p equal."""
        U_n = float(bsk.U_q(0.08, 0.08, 'n'))
        self.assertAlmostEqual(U_n, float(bsk.U_q(0.08, 0.08, 'p')), places=10)
        self.assertLess(U_n, -50)
        self.assertGreater(U_n, -80)


class TestGradientMeanField(unittest.TestCase):
    """U_LaplaceRho_q must be the functional derivative of int epsilon_delta_rho_np over rho_q,
    and equal U_LaplaceRho_n / U_LaplaceRho_p of W-BSK (hpc-engine/bsk/bsk_functional.c)."""

    # rho_n, rho_p, (grad rho_n)^2, (grad rho_p)^2, (grad rho)^2, lap rho_n, lap rho_p
    STATES = [(0.08, 0.07, 4e-4, 3e-4, 1.3e-3, -2e-3, -1.5e-3),
              (0.02, 0.005, 1e-4, 2e-5, 1.0e-4, 3e-3, 1e-3),
              (0.05, 0.0, 6e-4, 0.0, 6e-4, -4e-3, 0.0),
              (0.003, 0.001, 5e-5, 1e-5, 2e-5, 2e-3, -1e-3)]
    # U_LaplaceRho_n, U_LaplaceRho_p, e_LaplaceRho_np of W-BSK (BSk31, NUCLEAR_MATTER), computed
    # with bsk_functional.c through ctypes (2026-11-27); DENSEPSILON 1e-8 there, 1e-12 here.
    WBSK = [(0.4464958894, 0.4738608131, 6.960472192043e-02),
            (-0.4883970953, -0.5309509298, 3.310879766109e-03),
            (0.3496947845, 1.0351161976, 3.340271511628e-02),
            (0.0292249519, -0.1399124218, -3.576142271957e-03)]

    def check_wbsk(self, mod):
        for (rn, rp, gn, gp, g, ln, lp), (u_n, u_p, e) in zip(self.STATES, self.WBSK):
            with self.subTest(module=mod.__name__, rho_n=rn, rho_p=rp):
                self.assertAlmostEqual(float(mod.U_LaplaceRho_q(rn, rp, 'n', gn, gp, g, ln, lp)), u_n, delta=1e-5)
                self.assertAlmostEqual(float(mod.U_LaplaceRho_q(rn, rp, 'p', gn, gp, g, ln, lp)), u_p, delta=1e-5)
                self.assertAlmostEqual(float(mod.epsilon_delta_rho_np(rn, rp, gn, gp, g)), e, delta=1e-8)

    def check_functional_derivative(self, mod):
        """d/ds int e(rho + s h) dx = int U h dx on a periodic 1D grid (spectral derivatives)."""
        N, L = 256, 40.
        x = np.arange(N) * L / N
        k = 2 * np.pi * np.fft.fftfreq(N, L / N)
        D = lambda f: np.fft.ifft(1j * k * np.fft.fft(f)).real
        D2 = lambda f: np.fft.ifft(-k**2 * np.fft.fft(f)).real
        prof = lambda R, a, r0: r0 / (1 + np.exp((abs(x - L / 2) - R) / a))
        rho = {'n': prof(6, .6, .09) + 1e-3, 'p': prof(5.5, .5, .07) + 1e-4}
        h = np.exp(-((x - L / 2 - 4) / 1.5)**2)

        def E(rn, rp):
            gn, gp = D(rn), D(rp)
            return mod.epsilon_delta_rho_np(rn, rp, gn**2, gp**2, (gn + gp)**2).sum() * L / N

        gn, gp = D(rho['n']), D(rho['p'])
        for q in 'np':
            s = 1e-6
            plus, minus = dict(rho), dict(rho)
            plus[q], minus[q] = rho[q] + s * h, rho[q] - s * h
            num = (E(plus['n'], plus['p']) - E(minus['n'], minus['p'])) / (2 * s)
            U = mod.U_LaplaceRho_q(rho['n'], rho['p'], q, gn**2, gp**2, (gn + gp)**2,
                                   D2(rho['n']), D2(rho['p']))
            with self.subTest(module=mod.__name__, q=q):
                self.assertAlmostEqual((U * h).sum() * L / N, num, delta=1e-6 * abs(num))

    def test_bsk(self):
        self.check_wbsk(bsk)
        self.check_functional_derivative(bsk)

    def test_bskg(self):
        from libnest import bskg
        self.check_wbsk(bskg)
        self.check_functional_derivative(bskg)

    def test_U_q_with_gradients(self):
        rn, rp, gn, gp, g, ln, lp = self.STATES[0]
        for q in 'np':
            full = bsk.U_q(rn, rp, q, 0.1, 0.09, rho_grad_n_square=gn, rho_grad_p_square=gp,
                           rho_grad_square=g, rho_lap_n=ln, rho_lap_p=lp)
            self.assertAlmostEqual(float(full), float(bsk.U_q(rn, rp, q, 0.1, 0.09)
                                   + bsk.U_LaplaceRho_q(rn, rp, q, gn, gp, g, ln, lp)), places=12)
        with self.assertRaises(ValueError):
            bsk.U_q(rn, rp, 'n', rho_lap_n=ln, rho_lap_p=lp)

    def test_vacuum_is_zero(self):
        self.assertEqual(float(bsk.U_LaplaceRho_q(0., 0., 'n', 1e-4, 1e-4, 4e-4, 1e-3, 1e-3)), 0.)


if __name__ == '__main__':
    unittest.main()
