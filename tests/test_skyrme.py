#!/usr/bin/env python3
"""
Physics invariants of libnest.edf.SkyrmeFunctional (TODO_edf.md, phase E3).

These tests check relations that must hold for any parameter set — limits, derivatives
and symmetries — rather than re-implementing the formulas.
"""
import unittest

import numpy as np

from libnest import bsk, definitions
from libnest.edf import SkyrmeFunctional, SkyrmeParameters

# All Skyrme couplings zero: a free Fermi gas, M* = M.
FREE = SkyrmeFunctional(SkyrmeParameters(name="free gas", t0=0., t1=0., t2=0., t3=0.,
                                         x0=0., x1=0., t2x2=0., x3=0., alpha=1.))


class TestChemicalPotential(unittest.TestCase):

    def test_free_gas_mu_is_fermi_energy(self):
        """mu_q = hbar^2 kF^2 / 2M* reduces to the Fermi energy when M* = M."""
        rho_n = np.array([0.01, 0.05, 0.08])
        np.testing.assert_allclose(FREE.mu_q(rho_n, 0., "n"),
                                   definitions.eF_n(definitions.rho2kf(rho_n)), rtol=1e-9)

    def test_mu_scales_with_inverse_effective_mass(self):
        """mu_q / eF_q = M / M*_q."""
        rho_n, rho_p = 0.06, 0.02
        eF = definitions.eF_n(definitions.rho2kf(rho_n))
        self.assertAlmostEqual(bsk.mu_q(rho_n, rho_p, "n") / eF,
                               1. / bsk.effMn(rho_n, rho_p), places=9)

    def test_mixed_scalar_and_array_input(self):
        """An array for one density and a scalar for the other must work for both q."""
        rho_n = np.array([0.05, 0.08])
        mu_p = bsk.mu_q(rho_n, 0., "p")          # no protons: numerical zero
        np.testing.assert_array_equal(mu_p, [1e-12, 1e-12])
        mu_n = bsk.mu_q(rho_n, 0.01, "n")
        self.assertEqual(mu_n[1], bsk.mu_q(0.08, 0.01, "n"))
        self.assertEqual(bsk.I(rho_n, 0., "p").shape, (2,))

    def test_mu_magnitude_bsk31(self):
        """Neutron matter at 0.08 fm^-3: eF = 36.8 MeV and M*/M ~ 1, so mu = O(30-45) MeV."""
        self.assertTrue(25. < bsk.mu_q(0.08, 0., "n") < 45.)



class TestPairingStrength(unittest.TestCase):

    def test_gap_equation_normalization(self):
        """v_pi = -8 pi^2 / I_q * (hbar^2 / 2M*_q)^(3/2), with hbar^2/2M*_q = B_q."""
        for rho_n, rho_p, q in [(0.03, 0., "n"), (0.05, 0.02, "n"), (0.05, 0.02, "p")]:
            with self.subTest(rho_n=rho_n, rho_p=rho_p, q=q):
                ratio = (-bsk.v_pi(rho_n, rho_p, q) * bsk.I(rho_n, rho_p, q)
                         / (8. * np.pi**2 * bsk.B_q(rho_n, rho_p, q)**1.5))
                self.assertAlmostEqual(ratio, 1., places=9)

    def test_attractive_and_of_nuclear_size(self):
        """Pairing strength in neutron matter: attractive, a few hundred MeV fm^3."""
        v = bsk.v_pi(np.array([0.005, 0.01, 0.03, 0.05]), 0., "n")
        self.assertTrue(np.all(v < 0.))
        self.assertTrue(np.all(np.abs(v) < 2000.), v)



POINTS = [(0.08, 0.08), (0.06, 0.02), (0.02, 0.06), (0.1, 0.), (0.03, 0.01)]


def _d(f, x, h=1e-6):
    """Central finite difference df/dx."""
    return (f(x + h) - f(x - h)) / (2. * h)


class TestMeanFields(unittest.TestCase):
    """Mean fields are derivatives of the energy density."""

    def test_U_q_is_density_derivative_of_epsilon_rho(self):
        """U_q = d epsilon_rho / d rho_q (the t0 and t3 terms)."""
        for rho_n, rho_p in POINTS:
            with self.subTest(rho_n=rho_n, rho_p=rho_p):
                dn = _d(lambda x: bsk.epsilon_rho_np(x, rho_p), rho_n)
                dp = _d(lambda x: bsk.epsilon_rho_np(rho_n, x), rho_p)
                self.assertAlmostEqual(bsk.U_q(rho_n, rho_p, "n") / dn, 1., places=6)
                self.assertAlmostEqual(bsk.U_q(rho_n, rho_p, "p") / dp, 1., places=6)

    def test_B_q_is_tau_derivative_of_epsilon(self):
        """B_q = hbar^2/2m_q + d epsilon_tau / d tau_q."""
        p = bsk.BSK31.params
        tau_n, tau_p = 0.1, 0.05
        for rho_n, rho_p in POINTS:
            with self.subTest(rho_n=rho_n, rho_p=rho_p):
                dn = _d(lambda t: bsk.epsilon_tau_np(rho_n, rho_p, t, tau_p, 0., 0.), tau_n)
                dp = _d(lambda t: bsk.epsilon_tau_np(rho_n, rho_p, tau_n, t, 0., 0.), tau_p)
                self.assertAlmostEqual(bsk.B_q(rho_n, rho_p, "n"), p.hbar2m_n + dn, places=6)
                self.assertAlmostEqual(bsk.B_q(rho_n, rho_p, "p"), p.hbar2m_p + dp, places=6)


if __name__ == "__main__":
    unittest.main()
