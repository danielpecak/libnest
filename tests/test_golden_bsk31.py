#!/usr/bin/env python3
"""
Golden-master regression test for the BSk31 functional (TODO_edf.md, phase E0).

Re-evaluates every case of tests/golden/spec.py on the inputs stored in
tests/golden/bsk31.json and compares with the stored results, for array input and
point by point for scalar input (including which calls raise). It also checks that
the public API recorded in the file is still importable with compatible signatures.

A refactor must keep this test green without regenerating the file. For an intended
change of the numbers run ``python -m tests.golden.make_golden --write`` and explain
the change in the commit message.
"""
import importlib
import inspect
import json
import os
import unittest

import numpy as np

from tests.golden import spec

RTOL = 1e-12
ATOL = 1e-15  # far below NUMZERO = DENSEPSILON = 1e-12; only absorbs round-off around 0


class TestGoldenBSk31(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        path = os.path.join(os.path.dirname(spec.__file__), spec.GOLDEN_FILE)
        with open(path) as f:
            cls.golden = json.load(f)
        cls.inputs = {k: np.asarray(v, dtype=float) for k, v in cls.golden["inputs"].items()}

    def test_spec_and_golden_file_agree(self):
        """Every case in the spec has stored results, and vice versa."""
        self.assertEqual(set(spec.CASES), set(self.golden["cases"]))

    def test_results(self):
        for name, stored in self.golden["cases"].items():
            with self.subTest(case=name):
                case = {k: stored[k] for k in ("target", "args", "extra", "scalar")}
                new = spec.evaluate(case, self.inputs)
                self._compare_array(new, stored)
                if case["scalar"]:
                    self._compare_scalar(new, stored)

    def _compare_array(self, new, stored):
        self.assertEqual(new["mask"], stored["mask"], "mask of the array result")
        keep = ~np.asarray(stored["mask"], dtype=bool) if stored["mask"] is not None else ...
        np.testing.assert_allclose(np.asarray(new["array"])[keep],
                                   np.asarray(stored["array"])[keep],
                                   rtol=RTOL, atol=ATOL, equal_nan=True,
                                   err_msg="array input")

    def _compare_scalar(self, new, stored):
        self.assertEqual(new["scalar_error"], stored["scalar_error"],
                         "exceptions raised for scalar input")
        self.assertEqual(new["scalar_ndim"], stored["scalar_ndim"],
                         "dimension of the result for scalar input")
        as_array = lambda v: np.array([np.nan if x is None else x for x in v], dtype=float)
        np.testing.assert_allclose(as_array(new["scalar"]), as_array(stored["scalar"]),
                                   rtol=RTOL, atol=ATOL, equal_nan=True,
                                   err_msg="scalar input")

    def test_public_functions_compatible(self):
        """Recorded functions exist; old parameters are kept, new ones have defaults."""
        for module_name, functions in self.golden["api"]["functions"].items():
            module = importlib.import_module("libnest." + module_name)
            for name, old_params in functions.items():
                with self.subTest(function=module_name + "." + name):
                    self.assertTrue(hasattr(module, name), "missing")
                    params = inspect.signature(getattr(module, name)).parameters
                    self.assertEqual(list(params)[:len(old_params)], old_params)
                    for extra in list(params)[len(old_params):]:
                        self.assertIsNot(params[extra].default, inspect.Parameter.empty,
                                         "new parameter {!r} needs a default".format(extra))

    def test_public_constants_unchanged(self):
        for module_name, constants in self.golden["api"]["constants"].items():
            module = importlib.import_module("libnest." + module_name)
            for name, value in constants.items():
                with self.subTest(constant=module_name + "." + name):
                    self.assertEqual(float(getattr(module, name)), value)


if __name__ == "__main__":
    unittest.main()
