#!/usr/bin/env python3
"""
Tests for libnest.edf: parameter sets, the functional object and the registry.

The numbers of BSk31 are guarded by tests/test_golden_bsk31.py; here we test the
structure: lookup, immutability, independence of functionals, and that the
parameters are really used.
"""
import dataclasses
import unittest

import numpy as np

from libnest import bsk, definitions
from libnest.edf import (SkyrmeFunctional, available_functionals, get_functional,
                         register)


class TestRegistry(unittest.TestCase):

    def test_bsk31_lookup_is_case_insensitive(self):
        f = get_functional("BSk31")
        self.assertIs(f, bsk.BSK31)
        self.assertIs(get_functional("bsk31"), f)
        self.assertIs(get_functional("BSK31"), f)

    def test_available_functionals(self):
        self.assertIn("BSk31", available_functionals())

    def test_unknown_name_lists_the_available_ones(self):
        with self.assertRaises(KeyError) as ctx:
            get_functional("no-such-functional")
        self.assertIn("BSk31", str(ctx.exception))

    def test_register_is_idempotent_for_identical_sets(self):
        self.assertIs(register(bsk.BSK31.params, bsk.BSK31.pairing), bsk.BSK31)

    def test_register_refuses_a_different_set_under_a_taken_name(self):
        changed = dataclasses.replace(bsk.BSK31.params, t0=-2300.)
        with self.assertRaises(ValueError):
            register(changed, bsk.BSK31.pairing)


class TestFunctional(unittest.TestCase):

    def setUp(self):
        self.f = bsk.BSK31
        # t1 changes the effective mass, hence also mu_q and the pairing strength
        self.g = SkyrmeFunctional(
            dataclasses.replace(self.f.params, name="BSk31-modified", t1=700.),
            self.f.pairing)

    def test_parameters_are_immutable(self):
        with self.assertRaises(dataclasses.FrozenInstanceError):
            self.f.params.t0 = 0.

    def test_module_functions_are_the_bsk31_methods(self):
        self.assertEqual(bsk.energy_per_nucleon(0.08, 0.08),
                         self.f.energy_per_nucleon(0.08, 0.08))
        self.assertEqual(bsk.T0, self.f.params.t0)
        self.assertEqual(bsk.KAPPAN, self.f.pairing.kappa_n)

    def test_functionals_do_not_share_state(self):
        before = bsk.effMn(0.08, 0.0)
        modified = self.g.effMn(0.08, 0.0)
        self.assertNotEqual(modified, before)
        self.assertEqual(bsk.effMn(0.08, 0.0), before)

    def test_definitions_wrappers_take_a_functional(self):
        args = (0.05, 0.01, "n")
        self.assertEqual(definitions.mu_q(*args), self.f.mu_q(*args))
        self.assertEqual(definitions.mu_q(*args, functional="BSk31"), self.f.mu_q(*args))
        self.assertEqual(definitions.mu_q(*args, functional=self.g), self.g.mu_q(*args))
        self.assertNotEqual(self.g.mu_q(*args), self.f.mu_q(*args))

    def test_pairing_cutoff_is_a_parameter(self):
        h = SkyrmeFunctional(self.f.params, dataclasses.replace(self.f.pairing, cutoff=8.))
        self.assertNotEqual(h.I(0.02, 0.0, "n"), self.f.I(0.02, 0.0, "n"))

    def test_scalar_and_array_input(self):
        rho = np.array([0.02, 0.08])
        values = self.g.energy_per_nucleon(rho, rho)
        for i, r in enumerate(rho):
            self.assertAlmostEqual(values[i], self.g.energy_per_nucleon(r, r), places=12)

    def test_repr_and_name(self):
        self.assertEqual(self.f.name, "BSk31")
        self.assertEqual(repr(self.f), "SkyrmeFunctional('BSk31')")


if __name__ == "__main__":
    unittest.main()
