Inner Crust
===========
Tables and data concerning the inner crust, taken from multiple sources, for example :cite:p:`pecak2024WBSkMeff`.

.. contents:: Table of Contents
    :depth: 3
    :local:
    :backlinks: none

.. todo::
  Add information about Landau velocity

Bulk Neutron Properties
-----------------------
Quantities extracted from simulations for each density of the inner crust :math:`\bar\rho`: :math:`\rho_{Bn}` -- bulk density of neutrons, :math:`\Delta_n` -- pairing energy of neutrons, :math:`k_{\mathrm{F}}` -- wave vector calculated for bulk density of neutrons,  :math:`\epsilon_{\mathrm{F}}` -- Fermi energy, :math:`\epsilon_{\mathrm{F}}^*` -- Fermi energy calculated with respect to effective mass, :math:`N` -- number of neutrons, :math:`\xi` -- coherence length, :math:`R` -- radius of impurity, :math:`M_{\mathrm{eff}}` -- effective mass of impurity.

..  csv-table::
    :header: ":math:`\\bar\\rho` [fm :sup:`-3`]", ":math:`\\rho_{Bn}` [fm :sup:`-3`]", ":math:`\\Delta_n` [MeV]", ":math:`k_{\\mathrm{F}}` [fm :sup:`-1`]", ":math:`\\epsilon_{\\mathrm{F}}` [MeV]", ":math:`\\epsilon_{\\mathrm{F}}^*` [MeV]", ":math:`N`", ":math:`\\xi` [fm]", "R [fm]", ":math:`M_{\\mathrm{eff}}` [ :math:`m_n`]"
    :widths: 15, 12, 12, 15, 15, 15, 15, 15, 8, 15
    :file: data/inner_crust_bulk.csv


Effective Masses
----------------
Effective masses for dynamical, static, and hydro approaches for zirconium clusters in the inner crust. For the static case, one can think about :math:`M_{\mathrm{eff}}^s` semi-classically as the number of protons plus bound neutrons. However, due to the neutron medium, this is not the proper picture.

..  csv-table::
    :header: ":math:`\\bar\\rho` [fm :sup:`-3`]", ":math:`\\rho_{Bn}` [fm :sup:`-3`]",   ":math:`M_{\\mathrm{eff}}^d` [ :math:`m_n`]", ":math:`M_{\\mathrm{eff}}^s` [ :math:`m_n`]", ":math:`M_{\\mathrm{eff}}^h` [ :math:`m_n`]"
    :widths: 15, 15, 15, 15, 15
    :file: data/inner_crust_meff.csv



Collisions Initial State
------------------------
.. error::

  The system description is not accurate. Please correct it!

Two nuclei in an elongated box. The Coulomb interaction is present,
and the ultrarelativistic electrons are treated as a uniform gas.
Schematic:

.. code::

      +-----------+
      | o |   | o |
      +-----------+

The nuclei are two zirconium atoms, so in total there are :math:`Z=80` protons. The units of energy :math:`E_{\mathrm{tot}}`, pairing :math:`\Delta`, and chemical potential :math:`\mu` are in MeV. The densities are measured in fm :sup:`-3`.

..  csv-table::
    :header: ":math:`\\bar\\rho`", ":math:`\\rho_{Bn}`", ":math:`\\mu_p`", :math:`N` , ":math:`\\mu_n`", ":math:`\\Delta_n`", ":math:`\\Delta_p`", ":math:`E_{\\mathrm{tot}}`"
    :file: data/inner_crust_collisions.csv

