###########
Functionals
###########
libNeST implements the (extended) Skyrme energy density functional once, in
:class:`~libnest.edf.SkyrmeFunctional`. A parametrization is a set of numbers
(:class:`~libnest.edf.SkyrmeParameters` and :class:`~libnest.edf.PairingParameters`)
registered under a name. The page :doc:`edf` documents the API.

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none


Choosing a functional
=====================
.. code:: python

    from libnest.edf import available_functionals, get_functional

    available_functionals()                    # names of the registered functionals
    f = get_functional("BSk31")                # case-insensitive
    f.energy_per_nucleon(0.08, 0.08)           # E/A in symmetric matter [MeV]
    f.params.t0, f.pairing.scheme              # the numbers behind it

The functions of :mod:`libnest.bsk` are the methods of the BSk31 functional, so
``libnest.bsk.energy_per_nucleon`` and ``get_functional("BSk31").energy_per_nucleon`` are
the same. To vary a parametrization, derive a new parameter set with
:func:`dataclasses.replace` (see :doc:`edf`) or change the pairing part with
:meth:`~libnest.edf.SkyrmeFunctional.with_pairing` (see :doc:`pairing`).


Registered functionals
======================
The table is generated from the registry when the documentation is built.

.. include:: _generated/functionals.inc

The pairing parameters are explained on the page :doc:`pairing`. Please cite the reference
of the functional you use (see :doc:`help`).


Comparison
==========
.. figure:: _static/functionals.png
   :width: 100 %

   Energy per nucleon in symmetric matter and in pure neutron matter, and the neutron
   effective mass in neutron matter, for every registered functional.
