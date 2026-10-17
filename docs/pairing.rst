#########################
Pairing in uniform matter
#########################
In the BSk and BSkG functionals the strength of the :math:`{}^1S_0` pairing is not a free
parameter fitted to nuclei. It is fixed so that, in uniform matter, the functional
reproduces given *reference gaps* :math:`\Delta_q(\rho_n,\rho_p)`, :math:`q = n, p`
:cite:`chamel2008further,chamel2010effective`. Three ingredients define the pairing of a
functional; all three are selected in its :class:`~libnest.edf.PairingParameters` and
implemented in :mod:`libnest.edf.pairing`:

1. the reference gaps in pure neutron matter (NeuM) and in symmetric matter (SM),
2. an interpolation of these gaps to arbitrary asymmetry,
3. the cutoff :math:`\varepsilon_\Lambda` of the gap equation that turns the gaps into a
   pairing strength.

A gradient term, which vanishes in uniform matter, completes the pairing part of the
functional for inhomogeneous systems.

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none


Reference gaps
==============
The gaps of Cao, Lombardo and Schuck :cite:`cao2006screening`, calculated with medium
polarization and self-energy effects, are represented by the fits (gap model
``"cao2006"``)

.. math::

    \Delta_{\mathrm{NeuM}}(k_F) = \frac{3.37968\, k_F^2}{k_F^2+0.556092^2}
    \frac{(k_F-1.38236)^2}{(k_F-1.38236)^2+0.327517^2}, \qquad k_F < 1.38\ \mathrm{fm}^{-1},

    \Delta_{\mathrm{SM}}(k_F) = \frac{11.5586\, k_F^2}{k_F^2 + 0.489932^2}
    \frac{(k_F - 1.3142)^2}{(k_F - 1.3142)^2 + 0.906146^2}, \qquad k_F < 1.31\ \mathrm{fm}^{-1},

and zero above these Fermi momenta. :math:`k_F` is the Fermi wavevector of one species:
:math:`k_F = (3\pi^2\rho_n)^{1/3}` in neutron matter and :math:`(3\pi^2\rho/2)^{1/3}` in
symmetric matter of total density :math:`\rho`. The fits were made for BSk30–32
:cite:`chamel2016further` and are used by BSkG3 and BSkG4 as well (listed in Table 5 of
:cite:`grams2025skyrme`). They are the only gap model implemented so far.

.. image:: _static/pairing_vs_kf.png
  :width: 48 %
.. image:: _static/pairing_vs_rho.png
  :width: 48 %

See :func:`.neutron_pairing_field` and :func:`.symmetric_pairing_field`.


Asymmetric matter
=================
The many-body calculations give the gaps of NeuM and SM only; in between they are
interpolated. All three schemes in libNeST treat neutrons and protons alike,
:math:`\Delta_p(\rho_n,\rho_p) = \Delta_n(\rho_p,\rho_n)`, and are written here for species
:math:`q`, with :math:`\bar q` the other species and the asymmetry

.. math::

    \delta_q = \frac{\rho_q - \rho_{\bar q}}{\rho}, \qquad \rho = \rho_n + \rho_p .

``"chamel2009"`` — the BSk family since BSk17 :cite:`chamel2009pairing`:

.. math::

    \Delta_q = \Delta_{\mathrm{SM}}(\rho)\,(1 - |\delta_q|)
    + \delta_q \frac{\rho_q}{\rho} \Delta_{\mathrm{NeuM}}(\rho_q) .

For :math:`\delta_q` close to :math:`-1` (the minority species) this can turn negative;
libNeST sets the gap to zero there.

``"bskg3_eq3"`` — BSkG3 :cite:`grams2023skyrme`, positive by construction:

.. math::

    \Delta_q = \Delta_{\mathrm{SM}}(\rho)\,(1 - |\delta_q|) + |\delta_q|\, \Delta_{\mathrm{NeuM}}(\rho_q) .

``"bskg4_eq6"`` — BSkG4 :cite:`grams2025skyrme`, derived from the isospin splitting of the
effective mass:

.. math::

    \Delta_q = \Delta_{\mathrm{NeuM}}(\rho_q)
    \left[ \frac{\Delta_{\mathrm{SM}}(\rho)}{\Delta_{\mathrm{NeuM}}(\rho/2)} \right]^{1 - \delta_q} .

.. figure:: _static/pairing_schemes.png
   :width: 100 %

   Reference gaps of neutrons (top) and protons (bottom) for the three schemes, as a
   function of the total density at fixed proton fraction :math:`Y_p = \rho_p/\rho`.
   :math:`Y_p = 0.5` is symmetric matter (all schemes agree), :math:`Y_p = 0` pure neutron
   matter.

.. note::

   In the BSkG4 scheme the ratio :math:`\Delta_{\mathrm{SM}}(\rho)/\Delta_{\mathrm{NeuM}}(\rho/2)`
   vanishes where :math:`\Delta_{\mathrm{SM}}` has closed but :math:`\Delta_{\mathrm{NeuM}}(\rho/2)`
   has not (:math:`0.152 < \rho < 0.1775` fm :sup:`-3`); where both are closed it is taken as
   1, as in the C implementation of W-BSK. The proton gap of neutron-rich matter therefore
   drops to zero in that window and returns to :math:`\Delta_{\mathrm{NeuM}}(\rho_p)` above
   it (bottom right panel).

A functional uses the scheme named in its parameters; another one can be chosen without
changing the Skyrme part:

.. code:: python

    from libnest import bsk

    bsk.BSK31.pairing.scheme                       # 'chamel2009'
    g = bsk.BSK31.with_pairing(scheme="bskg4_eq6")
    g.neutron_ref_pairing_field(0.05, 0.01), bsk.neutron_ref_pairing_field(0.05, 0.01)

The functions :func:`.neutron_ref_pairing_field` and :func:`.proton_ref_pairing_field`
follow the functional's scheme; :func:`.neutron_ref_pairing_field_eq2`, ``_eq3`` and
``_eq6`` (and their proton counterparts) always use one scheme.


Pairing strength
================
Solving the gap equation of uniform matter for the strength that reproduces the reference
gap gives :cite:`chamel2010effective`

.. math::

    v^{\pi q}(\rho_n,\rho_p) = -\frac{8\pi^2}{I_q(\rho_n,\rho_p)}
    \left(\frac{\hbar^2}{2M^*_q(\rho_n,\rho_p)}\right)^{3/2},
    \qquad
    I_q = \sqrt{\mu_q} \left[ 2 \ln\frac{2\mu_q}{\Delta_q}
    + \Lambda\!\left(\frac{\varepsilon_\Lambda}{\mu_q}\right) \right],

with :math:`\mu_q = \hbar^2 k_{Fq}^2/2M^*_q` and the effective mass :math:`M^*_q` of the
functional, and the analytic approximation :math:`\Lambda(x)` of the remaining integral
(see also :cite:`pecak2021properties`). States are included up to :math:`\varepsilon_\Lambda`
above the chemical potential: 6.5 MeV for BSk31 :cite:`chamel2016further`, the parameter
:math:`E_{\mathrm{cut}}` of BSkG3 and BSkG4 (7.961 and 7.919 MeV). It is the ``cutoff`` of
:class:`~libnest.edf.PairingParameters`.

See :func:`.v_pi`, :func:`.I`, :func:`.Lambda` and :func:`.mu_q`.


Gradient term
=============
For inhomogeneous matter both families add a term in the gradient of the total density,
:math:`\rho = \rho_n + \rho_p`, to the pairing coupling :math:`g_q`, which enters the
pairing energy density :math:`\varepsilon_\pi = \frac{1}{4}\sum_q g_q \tilde\rho_q^2`
(:math:`\tilde\rho_q` is the anomalous density). The two families use different
conventions (``kappa_convention`` of :class:`~libnest.edf.PairingParameters`): the BSk
family adds the term to the strength (``"absolute"``, :math:`\kappa_q` in MeV fm :sup:`8`,
Eq. 9 of :cite:`chamel2016further`), the BSkG family scales the strength (``"relative"``,
:math:`\kappa_q` in fm :sup:`8`, Eq. 15 of :cite:`grams2023skyrme`):

.. math::

    g_q^{\mathrm{BSk}} = v^{\pi q} + \kappa_q (\nabla\rho)^2 ,
    \qquad
    g_q^{\mathrm{BSkG}} = v^{\pi q} \left[ 1 + \kappa_q (\nabla\rho)^2 \right] .

BSk31 has :math:`\kappa_n = -36630.4`, :math:`\kappa_p = -45207.2` MeV fm :sup:`8`
:cite:`chamel2016further`; BSkG3 and BSkG4 have :math:`\kappa_n = 123.20`,
:math:`\kappa_p = 129.07` fm :sup:`8` :cite:`grams2023skyrme,grams2025skyrme`. In both cases
the term makes pairing stronger at the surface of nuclei and clusters. In uniform matter
:math:`\nabla\rho = 0` and :math:`g_q = v^{\pi q}`.

See :func:`.g_q` and :func:`.epsilon_pi_np`.
