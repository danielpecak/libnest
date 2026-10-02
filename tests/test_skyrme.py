#!/usr/bin/env python3
"""
Physics invariants of libnest.edf.SkyrmeFunctional (TODO_edf.md, phase E3).

These tests check relations that must hold for any parameter set — limits, derivatives
and symmetries — rather than re-implementing the formulas.
"""
import dataclasses
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



class TestGradientTerms(unittest.TestCase):

    def _eps_np(self, rho_n, rho_p, g_n, g_p):
        """epsilon_np with uniform tau and no pairing (nu = 0)."""
        return bsk.epsilon_np(rho_n, rho_p, g_n, g_p, 0.1, 0.05, 0., 0., 0., 0.,
                              bsk.KAPPAN, bsk.KAPPAP)

    def test_epsilon_np_uses_the_total_density_gradient(self):
        """The gradient part of epsilon_np is epsilon_delta_rho_np with (grad rho)^2 =
        (grad rho_n + grad rho_p)^2, so it depends on the relative sign of the gradients."""
        rho_n, rho_p = 0.06, 0.02
        for g_n, g_p in [(0.01, 0.004), (0.01, -0.004)]:
            with self.subTest(g_n=g_n, g_p=g_p):
                grad_part = (self._eps_np(rho_n, rho_p, g_n, g_p)
                             - self._eps_np(rho_n, rho_p, 0., 0.))
                expected = bsk.epsilon_delta_rho_np(rho_n, rho_p, g_n**2, g_p**2, (g_n + g_p)**2)
                self.assertAlmostEqual(grad_part / expected, 1., places=9)


    def test_epsilon_delta_rho_is_symmetric_under_n_p_exchange(self):
        """The gradient energy is isospin symmetric: swapping n and p changes nothing."""
        for rho_n, rho_p, g_n, g_p in [(0.06, 0.02, 0.01, 0.004), (0.05, 0.03, 0.02, -0.01)]:
            with self.subTest(rho_n=rho_n, rho_p=rho_p, g_n=g_n, g_p=g_p):
                g = (g_n + g_p)**2
                np_order = bsk.epsilon_delta_rho_np(rho_n, rho_p, g_n**2, g_p**2, g)
                pn_order = bsk.epsilon_delta_rho_np(rho_p, rho_n, g_p**2, g_n**2, g)
                self.assertAlmostEqual(np_order / pn_order, 1., places=12)



class TestKineticTerm(unittest.TestCase):
    """The kinetic energy uses the hbar^2/2m of the parameter set."""

    def setUp(self):
        params = dataclasses.replace(FREE.params, hbar2m_n=20.0, hbar2m_p=21.0)
        self.gas = SkyrmeFunctional(params)

    def test_free_gas_energy_per_nucleon(self):
        """Free Fermi gas: E/A = 3/5 hbar^2 kF^2 / 2m, averaged over the species."""
        rho = np.array([0.02, 0.08, 0.16])
        kF_neum = definitions.rho2kf(rho)
        np.testing.assert_allclose(self.gas.energy_per_nucleon(rho, 0.),
                                   0.6*20.0*kF_neum**2, rtol=1e-9)
        kF_snm = definitions.rho2kf(rho/2)
        np.testing.assert_allclose(self.gas.energy_per_nucleon(rho/2, rho/2),
                                   0.6*0.5*(20.0 + 21.0)*kF_snm**2, rtol=1e-9)

    def test_free_gas_energy_density(self):
        """Free Fermi gas without gradients and pairing: epsilon = sum_q hbar^2/2m_q tau_q."""
        eps = self.gas.epsilon_np(0.06, 0.02, 0., 0., 0.1, 0.05, 0., 0., 0., 0., 0., 0.)
        self.assertAlmostEqual(eps, 20.0*0.1 + 21.0*0.05, places=12)

    def test_general_and_neutron_matter_formulas_agree(self):
        """energy_per_nucleon(rho, 0) and the analytic energy_per_nucleon_n(rho)."""
        rho = np.array([0.01, 0.05, 0.1, 0.2])
        # rtol: the general formula adds DENSEPSILON = 1e-12 fm^-3 to rho (1e-10 relative
        # at 0.01 fm^-3); a bare-mass kinetic term instead of hbar2m would differ by ~1e-8
        np.testing.assert_allclose(bsk.energy_per_nucleon(rho, 0.),
                                   bsk.energy_per_nucleon_n(rho), rtol=1e-9)


if __name__ == "__main__":
    unittest.main()
