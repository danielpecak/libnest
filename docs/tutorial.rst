.. _tutorial:

Tutorial
========

First steps
-----------
Densities are in fm :sup:`-3`, energies in MeV, wavevectors in fm :sup:`-1`.
Functions accept both Python scalars and NumPy arrays.

.. code:: python

    import numpy as np
    from libnest import units, definitions, bsk

    # 1. Unit conversion: fm^-3 -> g/cm^3
    rho = 0.16
    print(f"{rho} fm^-3 = {units.fm3togcm3(rho):.3e} g/cm^3")

    # 2. Fermi wavevector and Fermi energy of neutrons.
    #    rho2kf takes the density of ONE nucleon species (both spin components).
    rho_n = 0.08
    kF = definitions.rho2kf(rho_n)
    print(f"kF = {kF:.3f} fm^-1, eF = {definitions.eF_n(kF):.2f} MeV")

    # 3. Energy per nucleon of symmetric matter (BSk31); arrays work as well.
    rho = np.linspace(0.04, 0.20, 5)               # total density [fm^-3]
    print(bsk.energy_per_nucleon(rho/2, rho/2))    # [MeV]

    # 4. Neutron 1S0 pairing gap in pure neutron matter
    print(f"Delta_n = {bsk.neutron_pairing_field(0.02):.3f} MeV")

Analyze your own data
---------------------
The script below loads a 2D neutron-density map, computes the Fermi wavevector
and the neutron pairing gap at every grid point, and saves a figure. It
generates a small sample dataset itself, so it runs as-is:

.. code:: console

    $ python examples/analyze_data.py

To use your own data, write files in the same format and point the
``LIBNEST_DATA`` environment variable at their directory.

.. literalinclude:: ../examples/analyze_data.py
   :language: python
