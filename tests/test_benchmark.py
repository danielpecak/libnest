#!/usr/bin/env python3
"""
Benchmark tests: libnest output vs published BSk31 nuclear-matter values.

These guard the physics, not just the code: the numbers below are external
reference values, so a regression that silently changes the functional will
fail here. Reference values for the Brussels-Montreal BSk31 functional follow
Goriely, Chamel & Pearson, Phys. Rev. C 93, 034337 (2016), plus standard
nuclear-matter phenomenology.
"""

import unittest
import numpy as np

from libnest import bsk, definitions


class TestBSk31Benchmarks(unittest.TestCase):
    """Symmetric nuclear matter (SNM) and neutron matter reference points."""

    @classmethod
    def setUpClass(cls):
        # E/A curve for symmetric matter (rho_n = rho_p = rho/2)
        cls.rho = np.linspace(0.10, 0.22, 4001)
        cls.EA = bsk.energy_per_nucleon(cls.rho / 2, cls.rho / 2)

    def test_saturation_density(self):
        """SNM saturates near the BSk31 value n0 = 0.1587 fm^-3."""
        rho0 = self.rho[np.argmin(self.EA)]
        self.assertAlmostEqual(rho0, 0.1587, delta=0.003)

    def test_saturation_energy(self):
        """Energy per nucleon at saturation ~ a_v = -16.05 MeV."""
        self.assertAlmostEqual(self.EA.min(), -16.05, delta=0.3)

    def test_snm_fermi_momentum(self):
        """SNM Fermi momentum at saturation ~ 1.33 fm^-1 (per nucleon species)."""
        self.assertAlmostEqual(float(definitions.rho2kf(0.08)), 1.33, delta=0.01)

    def test_neutron_matter_fermi_momentum(self):
        """Pure neutron matter at n0 = 0.16 fm^-3: kF ~ 1.68 fm^-1."""
        self.assertAlmostEqual(float(definitions.rho2kf(0.16)), 1.68, delta=0.01)

    def test_neutron_pairing_gap_scale(self):
        """Neutron 1S0 gap: O(1 MeV) peak in the low-density (kF < 1) regime."""
        rn = np.linspace(1e-4, 0.09, 20001)
        gap = bsk.neutron_pairing_field(rn)
        peak = gap.max()
        kF_at_peak = float(definitions.rho2kf(rn[np.argmax(gap)]))
        self.assertTrue(1.0 < peak < 2.5, f"unexpected peak gap {peak:.3f} MeV")
        self.assertTrue(0.6 < kF_at_peak < 1.1, f"unexpected kF {kF_at_peak:.3f}")


if __name__ == "__main__":
    unittest.main()
