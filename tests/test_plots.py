#!/usr/bin/env python3
"""
Smoke test for libnest.plots: every public plotting function runs and its figure
renders (mathtext labels are parsed only when a figure is drawn).
"""
import inspect
import unittest
from unittest import mock

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from libnest import plots

# Arguments for every public function of libnest.plots.
CALLS = {
    "plot_energy_per_nucleon": (0.5, 0.5),
    "plot_energy_per_nucleon_both": (),
    "plot_pairing_field_n": (1., 0.),
    "plot_pairing_field_p": (1., 0.5),
    "plot_effective_mass_n": (1., 1.),
    "plot_effective_mass_p": (1., 1.),
    "plot_B_q": (1., 1., "n"),
    "plot_U_q": (1., 1., "n"),
    "plot_isoscalarM": (1., 1.),
    "plot_isovectorM": (1., 1.),
    "plot_pressure_n": (0.2,),
    "plot_speed_of_sound_n": (0.2,),
    "plot_v_landau": (0.1,),
    "plot_v_critical": (0.1,),
    "plot_v_sf": (10.,),
    "plot_epsilon_rho_np": (0.5, 0.5),
    "plot_epsilon_tau_np": (0.5, 0.5, 0.1, 0.05, 0.01, 0.001),
    "plot_epsilon_delta_rho_np": (0.5, 0.5, 0.01, 0.005, 0.015),
    "plot_epsilon_np": (0.5, 0.5, 0.01, 0.005, 0.1, 0.05, 0.01, 0.001, 0.01, 0.01,
                        -36630.4, -45207.2),
}


def _draw_open_figures(*args, **kwargs):
    for number in plt.get_fignums():
        plt.figure(number).canvas.draw()


class TestPlots(unittest.TestCase):

    def test_every_public_function_is_covered(self):
        public = {n for n, f in inspect.getmembers(plots, inspect.isfunction)
                  if f.__module__ == plots.__name__ and not n.startswith("_")}
        self.assertEqual(public, set(CALLS))

    def test_functions_run_and_render(self):
        for name, args in CALLS.items():
            with self.subTest(function=name):
                try:
                    with mock.patch.object(plt, "show", _draw_open_figures):
                        getattr(plots, name)(*args)
                finally:
                    plt.close("all")


if __name__ == "__main__":
    unittest.main()
