#!/usr/bin/env python3
"""
Tests for libnest.edf.pairing: gap models, interpolation schemes, and how a
functional selects them.

The BSk31 numbers are guarded by tests/test_golden_bsk31.py; here we test the
limits every scheme must reproduce and the wiring of PairingParameters.
"""
import dataclasses
import unittest

import numpy as np

from libnest import bsk
from libnest.definitions import rho2kf
from libnest.edf import PairingParameters, pairing
from libnest.units import NUMZERO

GAPS = pairing.GAP_MODELS["cao2006"]
# densities where both reference gaps are open (kF(rho/2) < 1.31 and kF(rho) < 1.38)
RHO = np.array([1e-3, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08])


class TestGapModel(unittest.TestCase):

    def test_gaps_close_above_their_cutoffs(self):
        for kF_max, gap in ((1.38, pairing.cao2006_neutron_matter),
                            (1.31, pairing.cao2006_symmetric_matter)):
            with self.subTest(gap=gap.__name__):
                self.assertGreater(gap(kF_max - 0.01), NUMZERO)
                self.assertEqual(gap(kF_max + 0.01), NUMZERO)
                np.testing.assert_array_equal(
                    gap(np.array([kF_max + 0.01, kF_max + 0.5])), [NUMZERO, NUMZERO])

    def test_density_forms_use_the_per_species_fermi_momentum(self):
        rho = 0.05
        self.assertEqual(GAPS.neutron_matter(rho), GAPS.neutron_matter_kf(rho2kf(rho)))
        self.assertAlmostEqual(GAPS.symmetric_matter(0.3*rho, 0.7*rho),
                               GAPS.symmetric_matter_kf(rho2kf(rho/2)), places=12)


class TestSchemes(unittest.TestCase):
    """Limits shared by every registered scheme."""

    def test_symmetric_matter_limit(self):
        for name, scheme in pairing.SCHEMES.items():
            with self.subTest(scheme=name):
                np.testing.assert_allclose(scheme(GAPS, RHO/2, RHO/2),
                                           GAPS.symmetric_matter(RHO/2, RHO/2), rtol=1e-9)

    def test_pure_matter_limit(self):
        for name, scheme in pairing.SCHEMES.items():
            with self.subTest(scheme=name):
                np.testing.assert_allclose(scheme(GAPS, RHO, 0.*RHO),
                                           GAPS.neutron_matter(RHO), rtol=1e-8)

    def test_non_negative_at_any_asymmetry(self):
        rho, x = np.meshgrid(np.linspace(1e-4, 0.3, 61), np.linspace(0., 1., 101))
        for name, scheme in pairing.SCHEMES.items():
            with self.subTest(scheme=name):
                delta = scheme(GAPS, x*rho, (1.-x)*rho)
                self.assertTrue(np.all(np.isfinite(delta)))
                self.assertTrue(np.all(delta >= NUMZERO))

    def test_scalar_and_array_input(self):
        rho_q, rho_qbar = np.array([0.03, 0.05]), np.array([0.01, 0.0])
        for name, scheme in pairing.SCHEMES.items():
            with self.subTest(scheme=name):
                values = scheme(GAPS, rho_q, rho_qbar)
                for i in range(len(rho_q)):
                    value = scheme(GAPS, float(rho_q[i]), float(rho_qbar[i]))
                    self.assertIsInstance(value, float)
                    self.assertEqual(value, values[i])


class TestFunctionalPairing(unittest.TestCase):

    def setUp(self):
        self.f = bsk.BSK31
        self.rho_n = np.array([0.01, 0.03, 0.05, 0.07])
        self.rho_p = np.array([0.005, 0.01, 0.002, 0.0])

    def test_bsk31_uses_the_bsk_scheme(self):
        self.assertEqual(self.f.pairing.scheme, "chamel2009")
        self.assertEqual(self.f.pairing.gap_model, "cao2006")

    def test_with_pairing_selects_the_scheme(self):
        for scheme, suffix in (("chamel2009", "eq2"), ("bskg3_eq3", "eq3"),
                               ("bskg4_eq6", "eq6")):
            g = self.f.with_pairing(scheme=scheme)
            for q in ("neutron", "proton"):
                with self.subTest(scheme=scheme, q=q):
                    np.testing.assert_array_equal(
                        getattr(g, q + "_ref_pairing_field")(self.rho_n, self.rho_p),
                        getattr(self.f, q + "_ref_pairing_field_" + suffix)(self.rho_n,
                                                                            self.rho_p))

    def test_with_pairing_leaves_the_original_unchanged(self):
        g = self.f.with_pairing(scheme="bskg4_eq6", cutoff=7.9)
        self.assertIs(g.params, self.f.params)
        self.assertEqual(self.f.pairing.scheme, "chamel2009")
        self.assertEqual(self.f.pairing.cutoff, 6.5)
        self.assertEqual(g.pairing.cutoff, 7.9)

    def test_pairing_strength_follows_the_scheme(self):
        g = self.f.with_pairing(scheme="bskg4_eq6")
        # asymmetric matter: different reference gaps, hence different strengths
        self.assertNotAlmostEqual(g.v_pi(0.03, 0.01, "n"), self.f.v_pi(0.03, 0.01, "n"))
        # symmetric matter: every scheme gives Delta_SM
        self.assertAlmostEqual(g.v_pi(0.02, 0.02, "p"), self.f.v_pi(0.02, 0.02, "p"),
                               places=6)

    def test_unknown_names_are_rejected(self):
        with self.assertRaises(ValueError):
            PairingParameters(scheme="no-such-scheme")
        with self.assertRaises(ValueError):
            PairingParameters(gap_model="no-such-model")
        with self.assertRaises(ValueError):
            self.f.with_pairing(scheme="Chamel2009")

    def test_parameters_stay_immutable(self):
        with self.assertRaises(dataclasses.FrozenInstanceError):
            self.f.pairing.scheme = "bskg4_eq6"


class TestKappaConvention(unittest.TestCase):
    """g_q = v_pi + kappa (grad rho)^2 (BSk, "absolute") or v_pi [1 + kappa (grad rho)^2]
    (BSkG, "relative")."""

    rho = (0.05, 0.02)
    grad2 = 1e-4        # (grad rho)^2 [fm^-8]

    def setUp(self):
        self.absolute = bsk.BSK31
        # BSkG3/BSkG4 gradient pairing (kappa in fm^8) on the Skyrme part of BSk31
        self.relative = bsk.BSK31.with_pairing(kappa_convention="relative",
                                               kappa_n=123.20, kappa_p=129.07)

    def test_bsk31_is_absolute(self):
        self.assertEqual(self.absolute.pairing.kappa_convention, "absolute")

    def test_uniform_matter_coupling_is_the_pairing_strength(self):
        for f in (self.absolute, self.relative):
            for q in "np":
                with self.subTest(convention=f.pairing.kappa_convention, q=q):
                    self.assertEqual(f.g_q(*self.rho, q), f.v_pi(*self.rho, q))

    def test_absolute_kappa_adds_to_the_strength(self):
        for q, kappa in (("n", bsk.KAPPAN), ("p", bsk.KAPPAP)):
            with self.subTest(q=q):
                shift = self.absolute.g_q(*self.rho, q, self.grad2) - bsk.v_pi(*self.rho, q)
                self.assertAlmostEqual(shift / (kappa*self.grad2), 1., places=9)

    def test_relative_kappa_scales_the_strength(self):
        for q, kappa in (("n", 123.20), ("p", 129.07)):
            with self.subTest(q=q):
                ratio = (self.relative.g_q(*self.rho, q, self.grad2)
                         / self.relative.v_pi(*self.rho, q))
                self.assertAlmostEqual(ratio, 1. + kappa*self.grad2, places=12)

    def test_gradient_term_strengthens_pairing_with_the_published_kappas(self):
        """kappa < 0 (BSk31, absolute) and kappa > 0 (BSkG, relative) both make the
        coupling more attractive at the nuclear surface."""
        for f in (self.absolute, self.relative):
            for q in "np":
                with self.subTest(convention=f.pairing.kappa_convention, q=q):
                    self.assertLess(f.g_q(*self.rho, q), 0.)
                    self.assertLess(f.g_q(*self.rho, q, self.grad2), f.g_q(*self.rho, q))

    def test_epsilon_pi_uses_the_coupling_of_the_functional(self):
        rho_n, rho_p, g_n, g_p, nu_n, nu_p = 0.05, 0.02, 0.01, 0.004, 0.01, 0.004
        grad2 = (g_n + g_p)**2
        for f in (self.absolute, self.relative):
            with self.subTest(convention=f.pairing.kappa_convention):
                expected = 0.25*(f.g_q(rho_n, rho_p, "n", grad2)*nu_n**2
                                 + f.g_q(rho_n, rho_p, "p", grad2)*nu_p**2)
                self.assertAlmostEqual(
                    f.epsilon_pi_np(rho_n, rho_p, g_n, g_p, nu_n, nu_p) / expected, 1.,
                    places=12)
        # without kappa arguments, the functional's own are used
        self.assertEqual(bsk.epsilon_pi_np(rho_n, rho_p, g_n, g_p, nu_n, nu_p),
                         bsk.epsilon_pi_np(rho_n, rho_p, g_n, g_p, nu_n, nu_p,
                                           bsk.KAPPAN, bsk.KAPPAP))

    def test_scalar_and_array_input(self):
        rho_n, rho_p = np.array([0.03, 0.05]), np.array([0.01, 0.02])
        values = self.relative.g_q(rho_n, rho_p, "n", self.grad2)
        for i in range(len(rho_n)):
            self.assertAlmostEqual(
                values[i], self.relative.g_q(float(rho_n[i]), float(rho_p[i]), "n",
                                             self.grad2), places=12)

    def test_unknown_convention_is_rejected(self):
        with self.assertRaises(ValueError):
            PairingParameters(kappa_convention="multiplicative")


if __name__ == "__main__":
    unittest.main()
