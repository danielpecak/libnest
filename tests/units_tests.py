#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Info: who, where, when
# @author Daniel Pęcak <Daniel.Pecak@pw.edu.pl>
# Warsaw Technical University, Université Libre de Bruxelles
# On leave: Institute of Physics, Polish Academy of Sciences, Warsaw
# March 2022, Brussels
# =========== Description
# Script for converting the density units for nuclear matter
# in the ranges typical for neutron stars.
# =========== Usage example
# $ ./units_test.py


import unittest
import libnest.units

class TestUnits(unittest.TestCase):
    def test_fm3togcm3_zero(self):
        """Test zero density"""
        actual   = libnest.units.fm3togcm3(0.0)
        expected = 0.0
        self.assertEqual(actual, expected)

    def test_fm3togcm3_one(self):
        """Test 1 fm^{-3} density"""
        actual   = libnest.units.fm3togcm3(1.0)
        expected = 1.67377585e15  # More precise value
        self.assertAlmostEqual(actual, expected, places=5)


if __name__ == '__main__':
    unittest.main()
