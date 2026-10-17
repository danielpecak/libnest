"""Shared helpers for the documentation figure scripts (plot_*.py).

Each script is run by docs/source/Makefile as ``python plots/plot_<name>.py <output.png>``.
"""
import sys

import matplotlib
matplotlib.use("Agg")  # render to files, no display needed
import matplotlib.pyplot as plt

# Categorical colors, assigned in this fixed order (never cycled): checked for colour-vision
# deficiency and normal-vision separation of neighbours; the 3rd and 4th have low contrast on
# white, so every multi-series figure also needs a legend and distinct line styles.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
LINESTYLES = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (1, 3)), (0, (8, 3))]


def output_path():
    """Returns the output image path given as the only command-line argument."""
    if len(sys.argv) != 2:
        sys.exit("Usage: {} OUTPUT.png".format(sys.argv[0]))
    return sys.argv[1]


def savefig(path):
    """Saves the current figure to ``path`` and closes it."""
    plt.savefig(path, bbox_inches="tight")
    plt.close()
