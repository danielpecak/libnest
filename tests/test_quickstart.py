#!/usr/bin/env python3
"""
The Quick Start example exists once, in examples/quickstart.py: the documentation includes
that file, and the README must show the same code. Both must keep working.
"""
import contextlib
import io
import os
import re
import runpy
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
QUICKSTART = os.path.join(ROOT, "examples", "quickstart.py")


class TestQuickstart(unittest.TestCase):

    def test_readme_shows_the_quickstart_file(self):
        with open(os.path.join(ROOT, "README.md"), encoding="utf-8") as f:
            match = re.search(r"## 🚀 Quick Start\n\n```python\n(.*?)```", f.read(), re.S)
        self.assertIsNotNone(match, "README: no python block under '## 🚀 Quick Start'")
        with open(QUICKSTART, encoding="utf-8") as f:
            self.assertEqual(match.group(1), f.read(),
                             "README Quick Start and examples/quickstart.py differ")

    def test_quickstart_runs(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            runpy.run_path(QUICKSTART, run_name="__main__")
        self.assertIn("Energy per nucleon", output.getvalue())


if __name__ == "__main__":
    unittest.main()
