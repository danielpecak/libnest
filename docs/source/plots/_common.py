"""Shared helpers for the documentation figure scripts (plot_*.py).

Each script is run by docs/source/Makefile as ``python plots/plot_<name>.py <output.png>``.
"""
import sys

import matplotlib
matplotlib.use("Agg")  # render to files, no display needed
import matplotlib.pyplot as plt


def output_path():
    """Returns the output image path given as the only command-line argument."""
    if len(sys.argv) != 2:
        sys.exit("Usage: {} OUTPUT.png".format(sys.argv[0]))
    return sys.argv[1]


def savefig(path):
    """Saves the current figure to ``path`` and closes it."""
    plt.savefig(path, bbox_inches="tight")
    plt.close()
