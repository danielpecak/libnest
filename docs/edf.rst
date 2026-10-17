Module: EDF
===========
.. py:module:: libnest.edf

Energy density functionals: parameter sets, one implementation of the formulas,
and a registry to pick a functional by name.

.. code:: python

    from libnest.edf import get_functional, available_functionals

    print(available_functionals())
    f = get_functional("BSk31")                # case-insensitive
    f.energy_per_nucleon(0.08, 0.08)           # symmetric matter, E/A [MeV]
    f.params.t0, f.pairing.cutoff              # the parameters behind it

A parametrization is data: :class:`~libnest.edf.SkyrmeParameters` (Skyrme part)
and :class:`~libnest.edf.PairingParameters` (pairing part). The formulas are
implemented once, in :class:`~libnest.edf.SkyrmeFunctional`. The functions of
:mod:`libnest.bsk` are the methods of its BSk31 instance, so their documentation
applies to every functional.

To study a modified parametrization, derive a new parameter set:

.. code:: python

    from dataclasses import replace
    from libnest import bsk
    from libnest.edf import SkyrmeFunctional

    g = SkyrmeFunctional(replace(bsk.BSK31.params, name="BSk31-mod", t1=700.),
                         bsk.BSK31.pairing)
    g.effMn(0.08, 0.0), bsk.effMn(0.08, 0.0)

The pairing part works the same way; :meth:`~libnest.edf.SkyrmeFunctional.with_pairing`
returns a functional with the same Skyrme part and other pairing parameters, e.g. another
interpolation of the reference gaps to asymmetric matter:

.. code:: python

    h = bsk.BSK31.with_pairing(scheme="bskg4_eq6")
    h.neutron_ref_pairing_field(0.05, 0.01), bsk.neutron_ref_pairing_field(0.05, 0.01)

Registry
--------
.. autofunction:: libnest.edf.get_functional

.. autofunction:: libnest.edf.available_functionals

.. autofunction:: libnest.edf.register

Parameter sets
--------------
.. autoclass:: libnest.edf.SkyrmeParameters

.. autoclass:: libnest.edf.PairingParameters

Functional
----------
.. autoclass:: libnest.edf.SkyrmeFunctional
   :members: with_pairing

Pairing
-------
.. automodule:: libnest.edf.pairing
   :members: GapModel, cao2006_neutron_matter, cao2006_symmetric_matter,
             chamel2009, bskg3_eq3, bskg4_eq6

.. autodata:: libnest.edf.pairing.GAP_MODELS
   :no-value:

.. autodata:: libnest.edf.pairing.SCHEMES
   :no-value:
