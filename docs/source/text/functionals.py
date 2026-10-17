#!/usr/bin/env python3
"""Writes the table of all functionals registered in libnest.edf (an rst fragment).

Run by docs/source/Makefile as ``python text/functionals.py OUTPUT.inc``.
"""
import sys
from fractions import Fraction

from libnest.edf import available_functionals, get_functional


def number(x):
    return "{:.10g}".format(x)


def exponent(x):
    """1/5 instead of 0.2 for the density exponents."""
    fraction = Fraction(x).limit_denominator(100)
    return str(fraction) if abs(float(fraction) - x) < 1e-12 else number(x)


def doi(f):
    if not f.params.doi:
        return f.params.reference or "—"
    return "{} (`doi:{} <https://doi.org/{}>`__)".format(f.params.reference, f.params.doi,
                                                        f.params.doi)


FM3A = "MeV fm :sup:`3+3α`"
ROWS = [  # (label, unit, value of a functional)
    (":math:`t_0`", "MeV fm :sup:`3`", lambda f: number(f.params.t0)),
    (":math:`t_1`", "MeV fm :sup:`5`", lambda f: number(f.params.t1)),
    (":math:`t_2`", "MeV fm :sup:`5`", lambda f: number(f.params.t2)),
    (":math:`t_3`", FM3A, lambda f: number(f.params.t3)),
    (":math:`t_4`", "MeV fm :sup:`5+3β`", lambda f: number(f.params.t4)),
    (":math:`t_5`", "MeV fm :sup:`5+3γ`", lambda f: number(f.params.t5)),
    (":math:`x_0`", "", lambda f: number(f.params.x0)),
    (":math:`x_1`", "", lambda f: number(f.params.x1)),
    (":math:`t_2 x_2`", "MeV fm :sup:`5`", lambda f: number(f.params.t2x2)),
    (":math:`x_3`", "", lambda f: number(f.params.x3)),
    (":math:`x_4`", "", lambda f: number(f.params.x4)),
    (":math:`x_5`", "", lambda f: number(f.params.x5)),
    (":math:`\\alpha`", "", lambda f: exponent(f.params.alpha)),
    (":math:`\\beta`", "", lambda f: exponent(f.params.beta)),
    (":math:`\\gamma`", "", lambda f: exponent(f.params.gamma)),
    (":math:`\\hbar^2/2m_n`", "MeV fm :sup:`2`", lambda f: number(f.params.hbar2m_n)),
    (":math:`\\hbar^2/2m_p`", "MeV fm :sup:`2`", lambda f: number(f.params.hbar2m_p)),
    ("gap model", "", lambda f: "``{}``".format(f.pairing.gap_model)),
    ("interpolation scheme", "", lambda f: "``{}``".format(f.pairing.scheme)),
    (":math:`\\varepsilon_\\Lambda`", "MeV", lambda f: number(f.pairing.cutoff)),
    (":math:`\\kappa` convention", "", lambda f: "``{}``".format(f.pairing.kappa_convention)),
    (":math:`\\kappa_n`", "MeV fm :sup:`8` (absolute), fm :sup:`8` (relative)",
     lambda f: number(f.pairing.kappa_n)),
    (":math:`\\kappa_p`", "as :math:`\\kappa_n`", lambda f: number(f.pairing.kappa_p)),
    ("reference", "", doi),
]


def table(names):
    functionals = [get_functional(name) for name in names]
    lines = [".. list-table:: Parameters of the registered functionals",
             "   :header-rows: 1",
             "   :stub-columns: 1",
             "",
             "   * - Parameter",
             "     - Unit"]
    lines += ["     - {}".format(f.name) for f in functionals]
    for label, unit, value in ROWS:
        lines.append("   * - {}".format(label))
        lines.append("     - {}".format(unit))
        lines += ["     - {}".format(value(f)) for f in functionals]
    return "\n".join(line.rstrip() for line in lines) + "\n"


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: {} OUTPUT.inc".format(sys.argv[0]))
    with open(sys.argv[1], "w", encoding="utf-8") as f:
        f.write(table(available_functionals()))
