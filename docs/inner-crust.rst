Inner crust
===========
Tables and data considering inner crust, taken from :cite:p:`pecak2024WBSkMeff`.

Bulk neutron properties
-----------------------
Quantities extracted from simulations for each density of the inner crust :math:`\bar\rho`: :math:`\rho_{Bn}` -- bulk density of neutrons, :math:`\Delta_n` -- pairing energy of neutrons, :math:`k_{\mathrm{F}}` -- wave vector calculated for bulk density of neutrons,  :math:`\epsilon_{\mathrm{F}}` -- Fermi energy, :math:`\epsilon_{\mathrm{F}}^*` -- Fermi energy calculated with regard to effective mass, :math:`N` -- number of neutrons, :math:`\xi` -- coherence length, :math:`R` -- radius of impurity, :math:`M_{\mathrm{eff}}` -- effective mass of impurity.

..  csv-table::
    :header: ":math:`\\bar\\rho` [fm :sup:`-3`]", ":math:`\\rho_{Bn}` [fm :sup:`-3`]", ":math:`\\Delta_n` [MeV]", ":math:`k_{\\mathrm{F}}` [fm :sup:`-1`]", ":math:`\\epsilon_{\\mathrm{F}}`" [MeV], ":math:`\\epsilon_{\\mathrm{F}}^*` [\MeV]", ":math:`N`", ":math:`\\xi` [fm]", "R [fm]", ":math:`M_{\\mathrm{eff}}` [ :math:`m_n`]"
    :widths: 15, 12, 12, 15, 15, 15, 15, 15, 8, 15

    0.0023, 0.0016, 0.826, 0.363, 2.723, 2.694, 376.6, 5.79, 5.32, 150.75 :math:`\pm` 1.8
    0.0058, 0.0045, 1.236, 0.512, 5.436, 5.301, 934.4, 5.47, 5.21, 139.30 :math:`\pm` 1.0
    0.0104, 0.0084, 1.483, 0.628, 8.165, 7.840, 1696.8, 5.58, 5.64, 164.70 :math:`\pm` 3.1
    0.0148, 0.0120, 1.562, 0.707, 10.37, 9.839, 2385.3, 5.98, 5.79, 155.25 :math:`\pm` 2.9
    0.0187, 0.0152, 1.557, 0.766, 12.15, 11.421, 3032.5, 6.49, 6.21, 174.45 :math:`\pm` 4.1
    0.0237, 0.0193, 1.621, 0.829, 14.24, 13.248, 3789.3, 6.75, 5.28, 168.55 :math:`\pm` 6.3
    0.0267, 0.0217, 1.566, 0.863, 15.42, 14.268, 4258.0, 7.27, 5.55, 168.80 :math:`\pm` 7.1
    0.0300, 0.0244, 1.514, 0.898, 16.70, 15.368, 4847.1, 7.82, 7.46, 171.85 :math:`\pm` 8.3
    0.0338, 0.0276, 1.467, 0.935, 18.10, 16.558, 5430.9, 8.40, 7.27, 166.15 :math:`\pm` 8.3
    0.0428, 0.0351, 1.327, 1.013, 21.28, 19.260, 6925.0, 10.1, 8.71, 150.50 :math:`\pm` 8.8
    0.0510, 0.0422, 1.097, 1.077, 24.05, 21.618, 8070.6, 13.0, 9.03, 150.60 :math:`\pm` 8.8

Effective masses
----------------
Effective masses for dynamical, static, and hydro approaches for zirconium cluster in the inner crust. For the static case one can think about :math:`M_{\mathrm{eff}}^s` semi-classically like the number of protons plus bound neutrons. But due to the neutron medium it is not proper picture.

..  csv-table::
    :header: ":math:`\\bar\\rho` [fm :sup:`-3`]", ":math:`\\rho_{Bn}` [fm :sup:`-3`]",   ":math:`M_{\\mathrm{eff}}^d` [ :math:`m_n`]", ":math:`M_{\\mathrm{eff}}^s` [ :math:`m_n`]", ":math:`M_{\\mathrm{eff}}^h` [ :math:`m_n`]"
    :widths: 15, 15, 15, 15, 15

    0.0023, 0.0016, 150.75, 145.02, 55.98
    0.0058, 0.0045, 139.30, 144.50, 57.46
    0.0104, 0.0084, 164.70, 172.34, 62.29
    0.0148, 0.0120, 155.25, 172.56, 59.58
    0.0187, 0.0152, 174.45, 197.12, 63.93
    0.0237, 0.0193, 168.55, 190.26, 58.17
    0.0267, 0.0217, 168.80, 195.16, 56.31
    0.0300, 0.0244, 171.85, 215.65, 55.52
    0.0338, 0.0276, 166.15, 216.84, 53.57
    0.0428, 0.0351, 150.50, 229.78, 42.70
    0.0510, 0.0422, 150.60, 161.96, 34.32


Collisions initial state
------------------------
.. error::

  The system description is not true. Correct it!

Two nuclei in an elongated box. The Coulomb interaction is present,
the ultrarelativistic electrons are treated as un uniform gas.
Schematic:

.. code::

      +-----------+
      | o |   | o |
      +-----------+

The nuclei are two zirconium atoms, so in total there are :math:`Z=80` protons. The units of energy :math:`E_{\mathrm{tot}}`, pairing :math:`\Delta`, and chemical potential :math:`\mu` are in [MeV]. The densities are measured in [fm :sup:`-3`].

..  csv-table::
    :header: ":math:`\\bar\\rho`", ":math:`\\rho_{Bn}`", ":math:`\\mu_p`", :math:`N` , ":math:`\\mu_n`", ":math:`\\Delta_n`", ":math:`\\Delta_p`", ":math:`E_{\\mathrm{tot}}`"

    0.000, 0.000101,  -25.255,   246.725, 0.373642, 0.192, 0.0069, -5.420
    0.002, 0.001624,  -25.828,   623.797, 1.937040, 0.832, 0.0494, -1.796
    0.006, 0.004533,  -29.776,  1370.09,  3.555630, 1.256, 0.1061,  0.659
    0.010, 0.008369,  -37.631,  2404.02,  4.962623, 1.478, 0.1066,  2.230
    0.015, 0.011945,  -38.414,  3322.71,  5.862966, 1.579, 0.1193,  3.153
    0.019, 0.015194,  -42.382,  4203.90,  6.589333, 1.599, 0.1057,  3.825
    0.024, 0.019125,  -45.139,  5209.08,  7.239645, 1.639, 0.1057,  4.463
    0.027, 0.021536,  -45.138,  5836.39,  7.591119, 1.591, 0.1051,  4.792
    0.030, 0.024601,  -47.516,  6636.53,  8.032960, 1.538, 0.1031,  5.174
    0.034, 0.027579,  -49.270,  7409.12,  8.423588, 1.512, 0.1013,  5.508
    0.043, 0.035360,  -52.454,  9401.53,  9.212281, 1.354, 0.0907,  6.241
    0.051, 0.041282,  -55.162, 10923.8,   9.842774, 1.106, 0.0795,  6.712
