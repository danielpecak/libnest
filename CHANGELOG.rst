Changelog
=========
Notable changes of libNeST, newest first. Before version 1.0 a minor version may change the
public API; every such change is listed here.


Unreleased
----------
Several functionals with one implementation of the formulas (branch ``edf-refactor``).

Added
~~~~~
- ``libnest.edf``: the (extended) Skyrme functional implemented once, in
  ``SkyrmeFunctional``; parameter sets as immutable data (``SkyrmeParameters``,
  ``PairingParameters``); a registry (``get_functional``, ``available_functionals``,
  ``register``). ``libnest.bsk`` keeps its API: its functions are the methods of the BSk31
  instance ``bsk.BSK31``.
- ``libnest.edf.pairing``: the reference gaps (gap model ``"cao2006"``) and the interpolation
  schemes ``"chamel2009"`` (BSk), ``"bskg3_eq3"`` (BSkG3) and ``"bskg4_eq6"`` (BSkG4). A
  functional uses the scheme of its ``PairingParameters``; ``with_pairing(scheme=...)``
  returns a functional with another one.
- Gradient pairing in both conventions, ``kappa_convention="absolute"`` (BSk) and
  ``"relative"`` (BSkG), and ``g_q``, the local pairing coupling.
- ``definitions.mu_q``, ``xiBCS`` and ``E_minigap_rho_n`` take an optional ``functional=``.
- Documentation: pages on pairing and on the functionals, a figure gallery, this changelog,
  a generated BibTeX entry for libNeST.

Changed — results differ from 0.1.2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Physics fixes; each is covered by a test of a physical invariant.

- ``mu_q`` used the ratio M*/M as if it were a mass in MeV; it is now
  :math:`\hbar^2 k_F^2/2M^*_q`. This changes ``I``, ``v_pi``, ``epsilon_pi_np`` and
  ``epsilon_np``.
- ``v_pi`` uses :math:`(\hbar^2/2M^*_q)^{3/2}`.
- ``U_q``: the :math:`t_3` term is the derivative of ``epsilon_rho_np``.
- ``epsilon_np``: the square of the total density gradient is
  :math:`(\nabla\rho_n + \nabla\rho_p)^2`.
- ``epsilon_delta_rho_np``: the cross term is :math:`\nabla\rho_n\cdot\nabla\rho_p`.
- ``energy_per_nucleon`` and ``epsilon_np`` take :math:`\hbar^2/2m` from the parameter set
  (for BSk31 a relative change of at most 2e-7).
- ``epsilon_pi_np``: the :math:`\kappa` term uses the gradient of the total density
  (Eq. 9 of Phys. Rev. C 93, 034337).
- ``mu_q`` accepts an array for one density and a scalar for the other.
- ``kappa_n`` and ``kappa_p`` of ``epsilon_pi_np`` and ``epsilon_np`` are optional; the
  default is the functional's own.

Removed
~~~~~~~
- ``bsk.testMe`` (a debug helper; its formula is now a test).
- ``plots.plot_epsilon``, ``plot_epsilon_tau``, ``plot_epsilon_delta`` and ``epsilon_test``,
  which called functions that no longer exist.


0.1.2 — 2026-10-13
------------------
- Added the pairing interpolation schemes ``neutron/proton_ref_pairing_field_eq2``,
  ``_eq3`` and ``_eq6`` (BSkG3, BSkG4) with tests, and a template module ``bskg``.
- Fixed: ``symmetric_pairing_field`` takes the Fermi momentum of one species
  (:math:`\rho/2`); interpolated gaps no longer go negative; the formula of the proton
  reference gap; ``rho2tau`` (factor :math:`(3\pi^2)^{2/3}`).
- Documentation rebuilt: figures generated from the code, corrected docstrings, units and
  bibliography, strict Sphinx build in CI.


0.1.1 — 2026-07-14
------------------
- Tests against published BSk31 values (saturation point, Fermi momenta, gap scale) and a
  unified test suite; fixes in ``tools``.


0.1.0 — 2026-07-14
------------------
- First release on PyPI.
