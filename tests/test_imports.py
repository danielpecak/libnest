#!/usr/bin/env python3
"""
Import smoke tests.

Every module in the ``libnest`` package must import cleanly. This guards against
syntax/indentation errors and broken cross-module imports that would otherwise
only surface at runtime (a single broken low-level module such as ``bsk`` takes
the whole dependency stack down with it).
"""

import importlib
import pkgutil

import unittest

import libnest


# Discover every submodule of the libnest package so new modules are covered
# automatically without editing this test.
LIBNEST_MODULES = [
    f"libnest.{name}"
    for _, name, _ in pkgutil.iter_modules(libnest.__path__)
]


class TestImports(unittest.TestCase):
    """Each libnest submodule imports without raising."""

    def test_package_imports(self):
        self.assertTrue(hasattr(libnest, "__version__"))

    def test_all_submodules_import(self):
        self.assertGreater(len(LIBNEST_MODULES), 0, "no submodules discovered")
        for module_name in LIBNEST_MODULES:
            with self.subTest(module=module_name):
                importlib.import_module(module_name)


if __name__ == "__main__":
    unittest.main()
