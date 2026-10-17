#!/usr/bin/env python3
"""Writes the figure gallery (an rst fragment): every documentation figure, its caption
(the docstring of its script) and a download link to the script.

Run by docs/source/Makefile as ``python text/gallery.py OUTPUT.inc FIGURE.png ...``.
"""
import ast
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def entry(figure):
    stem = os.path.splitext(figure)[0]
    script = "plot_{}.py".format(stem)
    with open(os.path.join(HERE, "..", "plots", script), encoding="utf-8") as f:
        caption = ast.get_docstring(ast.parse(f.read())) or stem
    caption = " ".join(caption.split())
    return "\n".join([
        ".. figure:: /_static/{}".format(figure),
        "   :width: 70 %",
        "",
        "   {} Script: :download:`{} </source/plots/{}>`.".format(caption, script, script),
        "",
    ])


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("Usage: {} OUTPUT.inc FIGURE.png ...".format(sys.argv[0]))
    with open(sys.argv[1], "w", encoding="utf-8") as f:
        f.write("\n".join(entry(figure) for figure in sys.argv[2:]))
