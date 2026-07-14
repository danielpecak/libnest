#!/usr/bin/env python3
"""
Tests for libnest.units (physical constants and unit conversions).

Consolidates the former tests/units_tests.py (which pytest never discovered
because of its name) and the assertions implied by the old tests/tests.py
print-script (conversions must round-trip and accept scalars, lists, arrays).
"""

import unittest

import numpy as np

from libnest import units


class TestDensityConversion(unittest.TestCase):
    """fm^-3 <-> g/cm^3."""

    def test_fm3togcm3_zero(self):
        self.assertEqual(units.fm3togcm3(0.0), 0.0)

    def test_fm3togcm3_one(self):
        # 1 fm^-3 of nucleons ~ 1.67e15 g/cm^3
        self.assertAlmostEqual(units.fm3togcm3(1.0), 1.67377585e15, places=5)

    def test_density_roundtrip(self):
        for rho in (0.0, 0.08, 0.16, 0.5):
            self.assertAlmostEqual(units.gcm3tofm3(units.fm3togcm3(rho)), rho, places=10)


class TestTemperatureConversion(unittest.TestCase):
    """K <-> MeV."""

    def test_temperature_roundtrip(self):
        for t in (0.0, 1.0, 10.0):
            self.assertAlmostEqual(units.KtoMev(units.MeVtoK(t)), t, places=10)

    def test_mev_to_kelvin_scale(self):
        # 1 MeV ~ 1.16e10 K
        self.assertAlmostEqual(units.MeVtoK(1.0), 1.16e10, delta=1e8)


class TestConversionInputTypes(unittest.TestCase):
    """Conversions accept scalars, Python lists, and NumPy arrays."""

    def test_accepts_scalar_list_array(self):
        values = [0.0, 1.0, 2.0]
        for func in (units.KtoMev, units.MeVtoK, units.fm3togcm3, units.gcm3tofm3):
            with self.subTest(func=func.__name__):
                scalar = func(2.0)
                self.assertTrue(np.isscalar(scalar) or np.ndim(scalar) == 0)
                from_list = np.asarray(func(values))
                from_array = func(np.array(values))
                self.assertEqual(from_list.shape, (3,))
                np.testing.assert_allclose(from_list, from_array)


if __name__ == "__main__":
    unittest.main()
