#!/usr/bin/env python3
"""
Tests for libnest.myio (data loading).

Each test writes a tiny dataset in the documented format into a temporary
directory and exercises the readers against it — no data files are stored in
the repository.
"""

import os
import tempfile
import unittest
from pathlib import Path

from libnest import myio


def write_dataset(data_dir, prefix="sample", nx=8, ny=8, nz=1):
    """Write <prefix>_info.txt and <prefix>_density.txt (x y rho) fixtures."""
    data_dir = Path(data_dir)
    with open(data_dir / f"{prefix}_info.txt", "w") as f:
        f.write("# grid info\n")
        for k, v in [("NX", nx), ("NY", ny), ("NZ", nz),
                     ("DX", 1.0), ("DY", 1.0), ("DZ", 1.0)]:
            f.write(f"{k} = {v}\n")
    with open(data_dir / f"{prefix}_density.txt", "w") as f:
        f.write("# x y rho_n\n")
        for i in range(nx):
            for j in range(ny):
                f.write(f"{i:.1f} {j:.1f} {0.01 * (i + j):.6e}\n")
    return data_dir


class TestDimension(unittest.TestCase):
    """dimension() counts non-flat axes."""

    def test_2d(self):
        self.assertEqual(myio.dimension(8, 8, 1), 2)

    def test_1d(self):
        self.assertEqual(myio.dimension(8, 1, 1), 1)

    def test_3d(self):
        self.assertEqual(myio.dimension(8, 8, 8), 3)

    def test_0d(self):
        self.assertEqual(myio.dimension(1, 1, 1), 0)


class TestReaders(unittest.TestCase):
    """readDimTxt() and txt2df() against generated fixtures."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.data_dir = write_dataset(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def test_read_grid_info(self):
        dims = myio.readDimTxt(str(self.data_dir) + os.sep, "sample")
        self.assertEqual(dims, [8, 8, 1])

    def test_read_grid_info_generalizes_beyond_2d(self):
        # regression: a 3D grid must not hard-exit (old code only allowed 2D)
        write_dataset(self.data_dir, prefix="cube", nx=4, ny=4, nz=4)
        dims = myio.readDimTxt(str(self.data_dir) + os.sep, "cube")
        self.assertEqual(dims, [4, 4, 4])

    def test_read_grid_info_defaults_missing_dims_to_one(self):
        # an info file that lists only NX must default NY, NZ to 1 (not crash)
        with open(self.data_dir / "line_info.txt", "w") as f:
            f.write("NX = 5\n")
        dims = myio.readDimTxt(str(self.data_dir) + os.sep, "line")
        self.assertEqual(dims, [5, 1, 1])

    def test_txt2df_columns_and_shape(self):
        df = myio.txt2df(str(self.data_dir / "sample"), "density", ["x", "y", "rho"])
        self.assertEqual(list(df.columns), ["x", "y", "rho"])
        self.assertEqual(len(df), 64)  # 8 x 8 grid

    def test_txt2df_parses_values(self):
        df = myio.txt2df(str(self.data_dir / "sample"), "density", ["x", "y", "rho"])
        self.assertTrue((df["rho"] >= 0).all())


if __name__ == "__main__":
    unittest.main()
