#!/usr/bin/env python3
"""
Tests for libnest.tools (grid/field utilities).

These operate on gridded simulation data laid out as [time, x, y, z] (or a bare
spatial array for threeSlice), so the fixtures below are small numpy arrays with
hand-computable answers.
"""

import unittest

import numpy as np

from libnest import tools


class TestParticleN(unittest.TestCase):
    """particleN sums the density over the spatial axes per time step."""

    def test_counts_uniform_density(self):
        rho = np.full((1, 4, 4, 4), 0.1)          # 64 cells * 0.1
        np.testing.assert_allclose(tools.particleN(rho), [6.4])

    def test_per_timestep(self):
        rho = np.ones((3, 2, 2, 2))               # 8 per step, 3 steps
        np.testing.assert_allclose(tools.particleN(rho), [8.0, 8.0, 8.0])


class TestCenterOfMass(unittest.TestCase):
    """centerOfMass returns grid-index coordinates per time step."""

    def test_point_mass(self):
        rho = np.zeros((1, 4, 4, 4))
        rho[0, 1, 2, 3] = 5.0                      # all mass at (1, 2, 3)
        np.testing.assert_allclose(tools.centerOfMass(rho), [[1.0, 2.0, 3.0]])

    def test_uniform_is_geometric_center(self):
        rho = np.ones((1, 4, 4, 4))                # center of a 0..3 grid is 1.5
        np.testing.assert_allclose(tools.centerOfMass(rho), [[1.5, 1.5, 1.5]])


class TestThreeSlice(unittest.TestCase):
    """threeSlice cuts central slices; works for 1D, 2D, 3D (2D was crashing)."""

    def test_1d(self):
        a = np.arange(5)
        out = tools.threeSlice(a)
        np.testing.assert_array_equal(np.asarray(out[0]), a)

    def test_2d_returns_two_central_slices(self):
        a = np.arange(16).reshape(4, 4)
        out = tools.threeSlice(a)                   # regression: used to raise TypeError
        self.assertEqual(len(out), 2)
        np.testing.assert_array_equal(np.asarray(out[0]), a[:, 2])
        np.testing.assert_array_equal(np.asarray(out[1]), a[2, :])

    def test_3d_returns_three_central_slices(self):
        a = np.arange(27).reshape(3, 3, 3)
        out = tools.threeSlice(a)
        self.assertEqual(len(out), 3)
        np.testing.assert_array_equal(np.asarray(out[0]), a[:, 1, 1])
        np.testing.assert_array_equal(np.asarray(out[1]), a[1, :, 1])
        np.testing.assert_array_equal(np.asarray(out[2]), a[1, 1, :])


class TestCondensationEnergy(unittest.TestCase):
    """Smoke test: condensationEnergy evaluates to a finite value per time step."""

    def test_runs_and_finite(self):
        rho = np.full((1, 4, 4, 4), 0.05)
        delta = np.full((1, 4, 4, 4), 1.0)
        e = tools.condensationEnergy(rho, delta)
        self.assertEqual(e.shape, (1,))
        self.assertTrue(np.all(np.isfinite(e)))


if __name__ == "__main__":
    unittest.main()
