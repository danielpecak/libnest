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


if __name__ == "__main__":
    unittest.main()
