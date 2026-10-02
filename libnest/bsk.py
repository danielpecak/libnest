# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Daniel Pęcak
#
# Original paper with formulas:
#   https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804
# NOTE: Table I is outdated; use the following data:
# forBSk31
# According to Goriely, Chamel, Pearson PRC 93 034337 (2016)
#================================
"""
Module: BSk
===========
The Brussels-Montreal (BSk) family of Skyrme functionals: energy density,
effective masses, pairing fields, etc., for uniform matter.

The module-level functions below are the methods of the BSk31 functional
:data:`BSK31`. The formulas are implemented once, in
:class:`~libnest.edf.SkyrmeFunctional`; to use another parametrization, see
:func:`libnest.edf.get_functional`.


BSk31 Skyrme Parameters
-----------------------

The following table lists all Skyrme parameters used in the BSk31 functional:

.. list-table:: BSk31 Skyrme Parameters
   :header-rows: 1
   :widths: 15 15 40

   * - Parameter
     - Value
     - Unit
   * - t₀
     - -2302.01
     - MeV fm³
   * - t₁
     - 762.99
     - MeV fm⁵
   * - t₂
     - 0.0
     - MeV fm⁵
   * - t₃
     - 13797.83
     - MeV fm⁽³⁺³ᵅ⁾
   * - t₄
     - -500.0
     - MeV fm⁽⁵⁺³ᵝ⁾
   * - t₅
     - -40.0
     - MeV fm⁽⁵⁺³ᵞ⁾
   * - x₀
     - 0.676655
     - dimensionless
   * - x₁
     - 2.658109
     - dimensionless
   * - x₂t₂
     - -422.29
     - MeV fm⁵
   * - x₃
     - 0.83982
     - dimensionless
   * - x₄
     - 5.0
     - dimensionless
   * - x₅
     - -12.0
     - dimensionless
   * - α
     - 1/5
     - dimensionless
   * - β
     - 1/12
     - dimensionless
   * - γ
     - 1/4
     - dimensionless

Individual Parameters
---------------------

.. data:: T0   =-2302.01

    Skyrme parameter :math:`t_0` [MeV fm :sup:`3`]

.. data:: T1   =762.99

    Skyrme parameter :math:`t_1` [MeV fm :sup:`5`]

.. data:: T2   =0.0

    Skyrme parameter :math:`t_2` [MeV fm :sup:`5`]

.. data:: T3   =13797.83

    Skyrme parameter :math:`t_3` [MeV fm :sup:`(3+3*ALPHA)`]

.. data:: T4   =-500.

    Skyrme parameter :math:`t_4` [MeV fm :sup:`(5+3*BETA)`]

.. data:: T5   =-40.

    Skyrme parameter :math:`t_5` [MeV fm :sup:`(5+3*GAMMA)`]

.. data:: X0   =0.676655

    Skyrme parameter :math:`x_0` [1]

.. data:: X1   =2.658109

    Skyrme parameter :math:`x_1` [1]

.. data:: T2X2 =-422.29

    Skyrme parameter :math:`x_2t_2` [1][MeV fm :sup:`5`]

.. data:: X3   =0.83982

    Skyrme parameter :math:`x_3` [1]

.. data:: X4   =5.

    Skyrme parameter :math:`x_4` [1]

.. data:: X5   =-12.

    Skyrme parameter :math:`x_5` [1]

.. data:: ALPHA =(1./5.)

    [1]

.. data:: BETA  =(1./12.)

    [1]

.. data:: GAMMA =(1./4.)

    [1]

Note:
    these are not important at the moment

.. data:: YW   =2.

    [1]

.. data:: FNP  =1.00

    [1]

.. data:: FNM  =1.06

    [1]

.. data:: FPP  =1.00

    [1]

.. data:: FPM  =1.04

    [1]

.. data:: KAPPAN =-36630.4

    [MeV fm :sup:`8`]

.. data:: KAPPAP =-45207.2

    [MeV fm :sup:`8`]

List of functions
-----------------
"""

# Re-exported for backward compatibility: these names used to be importable from here.
from libnest.units import HBARC, DENSEPSILON, NUMZERO  # noqa: F401
from libnest.units import MN, MP, HBAR2M_n, HBAR2M_p  # noqa: F401
from libnest.definitions import rho2kf, rhoEta  # noqa: F401
from libnest.edf.parameters import SkyrmeParameters, PairingParameters
from libnest.edf.parametrizations import register
from libnest.edf.skyrme import Lambda

#: The BSk31 functional (:class:`~libnest.edf.skyrme.SkyrmeFunctional`); the
#: module-level functions of :mod:`libnest.bsk` are its methods.
BSK31 = register(
    SkyrmeParameters(
        name="BSk31",
        t0=-2302.01, t1=762.99, t2=0.0, t3=13797.83, t4=-500., t5=-40.,
        x0=0.676655, x1=2.658109, t2x2=-422.29, x3=0.83982, x4=5., x5=-12.,
        alpha=1./5., beta=1./12., gamma=1./4.,
        yw=2.,
        reference="S. Goriely, N. Chamel, J. M. Pearson, Phys. Rev. C 93, 034337 (2016)",
        doi="10.1103/PhysRevC.93.034337"),
    PairingParameters(cutoff=6.5, kappa_n=-36630.4, kappa_p=-45207.2,
                      f_np=1.00, f_nm=1.06, f_pp=1.00, f_pm=1.04))

# ================================
#   Module constants (BSk31), kept for backward compatibility
# ================================
T0    = BSK31.params.t0      # Skyrme parameter t0 [MeV fm^3]
T1    = BSK31.params.t1      # Skyrme parameter t1 [MeV fm^5]
T2    = BSK31.params.t2      # Skyrme parameter t2 [MeV fm^5]
T3    = BSK31.params.t3      # Skyrme parameter t3 [MeV fm^(3+3*ALPHA)]
T4    = BSK31.params.t4      # Skyrme parameter t4 [MeV fm^(5+3*BETA)]
T5    = BSK31.params.t5      # Skyrme parameter t5 [MeV fm^(5+3*GAMMA)]
X0    = BSK31.params.x0      # Skyrme parameter x0 [1]
X1    = BSK31.params.x1      # Skyrme parameter x1 [1]
T2X2  = BSK31.params.t2x2    # Skyrme parameter x2*t2 [MeV fm^5]
X3    = BSK31.params.x3      # Skyrme parameter x3 [1]
X4    = BSK31.params.x4      # Skyrme parameter x4 [1]
X5    = BSK31.params.x5      # Skyrme parameter x5 [1]
ALPHA = BSK31.params.alpha   # [1]
BETA  = BSK31.params.beta    # [1]
GAMMA = BSK31.params.gamma   # [1]
YW    = BSK31.params.yw      # [1]
FNP   = BSK31.pairing.f_np   # [1]
FNM   = BSK31.pairing.f_nm   # [1]
FPP   = BSK31.pairing.f_pp   # [1]
FPM   = BSK31.pairing.f_pm   # [1]
KAPPAN = BSK31.pairing.kappa_n  # [MeV fm^8]
KAPPAP = BSK31.pairing.kappa_p  # [MeV fm^8]

# ================================
#   Module-level API: the methods of BSK31
# ================================
neutron_pairing_field = BSK31.neutron_pairing_field
symmetric_pairing_field = BSK31.symmetric_pairing_field
neutron_ref_pairing_field = BSK31.neutron_ref_pairing_field
proton_ref_pairing_field = BSK31.proton_ref_pairing_field
neutron_ref_pairing_field_eq2 = BSK31.neutron_ref_pairing_field_eq2
proton_ref_pairing_field_eq2 = BSK31.proton_ref_pairing_field_eq2
neutron_ref_pairing_field_eq3 = BSK31.neutron_ref_pairing_field_eq3
proton_ref_pairing_field_eq3 = BSK31.proton_ref_pairing_field_eq3
neutron_ref_pairing_field_eq6 = BSK31.neutron_ref_pairing_field_eq6
proton_ref_pairing_field_eq6 = BSK31.proton_ref_pairing_field_eq6
energy_per_nucleon = BSK31.energy_per_nucleon
numerical_derivative_energy_per_nucleon_n = BSK31.numerical_derivative_energy_per_nucleon_n
numerical_second_derivative_energy_per_nucleon_n = BSK31.numerical_second_derivative_energy_per_nucleon_n
numerical_derivative_epsilon_n = BSK31.numerical_derivative_epsilon_n
numerical_pressure_n = BSK31.numerical_pressure_n
numerical_derivative_pressure_n = BSK31.numerical_derivative_pressure_n
numerical_speed_of_sound_n = BSK31.numerical_speed_of_sound_n
energy_per_nucleon_n = BSK31.energy_per_nucleon_n
derivative_energy_per_nucleon_n = BSK31.derivative_energy_per_nucleon_n
second_derivative_energy_per_nucleon_n = BSK31.second_derivative_energy_per_nucleon_n
pressure_n = BSK31.pressure_n
pressure_derivative_n = BSK31.pressure_derivative_n
epsilon_derivative_n = BSK31.epsilon_derivative_n
speed_of_sound_n = BSK31.speed_of_sound_n
isoscalarM = BSK31.isoscalarM
isovectorM = BSK31.isovectorM
effMn = BSK31.effMn
effMp = BSK31.effMp
U_q = BSK31.U_q
B_q = BSK31.B_q
C0_rho = BSK31.C0_rho
C1_rho = BSK31.C1_rho
C0_tau = BSK31.C0_tau
C1_tau = BSK31.C1_tau
epsilon_rho_np = BSK31.epsilon_rho_np
epsilon_tau_np = BSK31.epsilon_tau_np
epsilon_delta_rho_np = BSK31.epsilon_delta_rho_np
epsilon_np = BSK31.epsilon_np
epsilon_pi_np = BSK31.epsilon_pi_np
v_pi = BSK31.v_pi
I = BSK31.I

# used to be importable from here (it was imported from libnest.definitions)
mu_q = BSK31.mu_q

__all__ = (["BSK31", "T0", "T1", "T2", "T3", "T4", "T5", "X0", "X1", "T2X2", "X3", "X4", "X5",
            "ALPHA", "BETA", "GAMMA", "YW", "FNP", "FNM", "FPP", "FPM", "KAPPAN", "KAPPAP"]
           + [
    "neutron_pairing_field",
    "symmetric_pairing_field",
    "neutron_ref_pairing_field",
    "proton_ref_pairing_field",
    "neutron_ref_pairing_field_eq2",
    "proton_ref_pairing_field_eq2",
    "neutron_ref_pairing_field_eq3",
    "proton_ref_pairing_field_eq3",
    "neutron_ref_pairing_field_eq6",
    "proton_ref_pairing_field_eq6",
    "energy_per_nucleon",
    "numerical_derivative_energy_per_nucleon_n",
    "numerical_second_derivative_energy_per_nucleon_n",
    "numerical_derivative_epsilon_n",
    "numerical_pressure_n",
    "numerical_derivative_pressure_n",
    "numerical_speed_of_sound_n",
    "energy_per_nucleon_n",
    "derivative_energy_per_nucleon_n",
    "second_derivative_energy_per_nucleon_n",
    "pressure_n",
    "pressure_derivative_n",
    "epsilon_derivative_n",
    "speed_of_sound_n",
    "isoscalarM",
    "isovectorM",
    "effMn",
    "effMp",
    "U_q",
    "B_q",
    "C0_rho",
    "C1_rho",
    "C0_tau",
    "C1_tau",
    "epsilon_rho_np",
    "epsilon_tau_np",
    "epsilon_delta_rho_np",
    "epsilon_np",
    "epsilon_pi_np",
    "v_pi",
    "I",
    "Lambda",
])
