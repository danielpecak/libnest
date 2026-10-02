# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Daniel Pęcak
"""
One implementation of the (extended) Skyrme functional for uniform matter.

:class:`SkyrmeFunctional` holds the formulas; the numbers come from its
:class:`~libnest.edf.parameters.SkyrmeParameters` and
:class:`~libnest.edf.parameters.PairingParameters`. The module-level functions of
:mod:`libnest.bsk` are the methods of the BSk31 instance.

Formulas and conventions as in :mod:`libnest.bsk`: densities in fm :sup:`-3` (one
nucleon species, both spin components), energies in MeV; scalars and NumPy arrays
are both accepted.
"""
import sys
import numpy as np
from libnest.units import HBARC, DENSEPSILON, NUMZERO
from libnest.units import MN, MP
from libnest.definitions import rho2kf, rhoEta
from libnest.edf.parameters import PairingParameters


def _clamp_nonnegative(pairing):
    """Floors an interpolated pairing field at (numerical) zero.

    The SM/NeuM interpolations carry a minus sign in one of the two terms, so
    for sufficiently asymmetric matter the result can dip below zero.

    Args:
        pairing (float or np.ndarray): raw interpolated pairing field [MeV]

    Returns:
        float or np.ndarray: the same, with negative values replaced by NUMZERO
    """
    pairing = np.maximum(np.asarray(pairing, dtype=float), NUMZERO)
    if pairing.shape == ():
        return float(pairing)
    return pairing


def Lambda(x):
    # Equation 16 from Phys Rev C 104
    """
    Function :math:`\\Lambda` used in :func:`.I`.
    Check :cite:`pecak2021properties` and earlier papers of Chamel.

    .. math::

        \\Lambda(x)
        = \\ln(16x) + 2\\sqrt{1+x} - 2 \\ln\\left({1+\\sqrt{1+x}}\\right) - 4.

    Args:
        x (float): variable

    Returns:
        float: `\\Lambda`

    See also:
        :func:`.I`
    """
    x = np.asarray(x)
    # Use np.where for vectorized, robust computation
    result = np.where(
        x <= 0,
        NUMZERO,
        np.log(16. * x) + 2. * np.sqrt(1. + x) - 2. * np.log(1. + np.sqrt(1. + x)) - 4.
    )
    # Return scalar if input was scalar
    if result.shape == ():
        return float(result)
    return result


class SkyrmeFunctional:
    """
    A Skyrme-type energy density functional: parameters plus the formulas.

    The methods have the names and signatures of the functions in
    :mod:`libnest.bsk` (which are the methods of the BSk31 instance).

    Args:
        params (SkyrmeParameters): Skyrme part
        pairing (PairingParameters): pairing part; defaults to ``PairingParameters()``

    Example:
        >>> from libnest.edf import get_functional
        >>> f = get_functional("BSk31")
        >>> e = f.energy_per_nucleon(0.08, 0.08)
    """

    def __init__(self, params, pairing=None):
        self.params = params
        self.pairing = PairingParameters() if pairing is None else pairing

    @property
    def name(self):
        """str: name of the parametrization"""
        return self.params.name

    def __repr__(self):
        return "SkyrmeFunctional({!r})".format(self.name)

    # ================================
    #       Pairing fields
    # ================================
    def neutron_pairing_field(self, rho_n):
        #    Formula (5.12) from NeST.pdf
        r"""
        Returns the pairing field for uniform pure neutron nuclear matter. For kF larger
        than 1.38 fm :sup:`-1` it returns (numerical) zero.

        .. math::

    	   \Delta_{\mathrm{NeuM}}(k_F) = \frac{3.37968 k_F^2}{k_F^2+0.556092^2} \frac{(k_F-1.38236)^2}{(k_F-1.38236)^2+0.327517^2},

        Args:
            rho_n (float): neutron density :math:`\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: pairing field for neutron matter :math:`\Delta_{\mathrm{NeuM}}` [MeV]

        See also:
            :func:`.symmetric_pairing_field`
            :func:`.neutron_ref_pairing_field`
            :func:`.proton_ref_pairing_field`
        """
        kF = np.asarray(rho2kf(rho_n), dtype=float)
        delta = 3.37968*(kF**2)*((kF-1.38236)**2)/(((kF**2)+(0.556092**2))*
                                                  ((kF-1.38236)**2+(0.327517**2)))
        if delta.shape == ():
            if kF > 1.38:
                return float(NUMZERO)
            else:
                return float(delta)
        delta[np.where(kF>1.38)] = np.float64(1.0*NUMZERO)
        return delta

    def symmetric_pairing_field(self, rho_n, rho_p):
        #   Formula (5.11) from NeST.pdf
        r"""
        Returns the pairing field for uniform symmetric matter. For kF larger than
        1.31 fm :sup:`-1` it returns (numerical) zero.

        .. math::

    	   \Delta_{\mathrm{SM}}(k_F) =  \frac{11.5586 k_F^2}{k_F^2 + 0.489932^2}\frac{(k_F - 1.3142)^2}{(k_F - 1.3142)^2 + 0.906146^2}.

        Args:
            rho_n (float): neutron density :math:`\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: pairing field for symmetric matter :math:`\Delta_{\mathrm{SM}}` [MeV]


        See also:
            :func:`.neutron_pairing_field`
            :func:`.neutron_ref_pairing_field`
            :func:`.proton_ref_pairing_field`
        """
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        #   The parametrization is a function of the PER-SPECIES Fermi momentum: in
        #   symmetric matter rho_n = rho_p = rho/2, hence the 0.5 factor (same
        #   convention as formula (A14) used in energy_per_nucleon). 
        #    https://doi.org/10.1103/PhysRevC.80.065804
        #    Chamel, N., Goriely, S., & Pearson, J. M. (2009). PRC 80, 065804 (2009).
        kF = np.asarray(rho2kf(0.5*(rho_n+rho_p)), dtype=float)
        delta = 11.5586*(kF**2)*((kF-1.3142)**2)/(((kF**2)+(0.489932**2))*
                                                 (((kF-1.3142)**2)+(0.906146**2)))
        if delta.shape == ():
            if kF > 1.31:
                return float(NUMZERO)
            else:
                return float(delta)
        delta[np.where(kF>1.31)] = NUMZERO
        return delta


    def neutron_ref_pairing_field(self, rho_n, rho_p):
        #   Formula (5.10) from NeST.pdf
        r"""
        Returns the reference pairing field for neutrons in uniform matter.
        This is an extrapolation between :math:`\Delta_{\mathrm{SM}}` and
        :math:`\Delta_{\mathrm{NeuM}}`. In limits :math:`\eta \rightarrow 0` reproduces
        symmetric matter and :math:`\eta \rightarrow 1`, the neutron matter.

        .. math::

            \Delta_n(\rho_n,\rho_p) =
            \Delta_{\mathrm{SM}}(\rho_n+\rho_p)
            \left( 1 - |\eta| \right)
            + \Delta_{\mathrm{NeuM}}(\rho_n) \eta \frac{\rho_n}{\rho_n+\rho_p}


        Args:
            rho_n (float): neutron density :math:`\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: pairing field for neutrons :math:`\Delta_n` [MeV]

        See also:
            :func:`.neutron_pairing_field`
            :func:`.symmetric_pairing_field`
            :func:`.proton_ref_pairing_field`
        """
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        rho, eta = rhoEta(rho_n, rho_p)
        rho = np.asarray(rho + DENSEPSILON, dtype=float)
        return _clamp_nonnegative(self.symmetric_pairing_field(rho_n, rho_p)*(1-np.abs(eta/rho))
            +self.neutron_pairing_field(rho_n)*rho_n/rho*eta/rho)

    def proton_ref_pairing_field(self, rho_n, rho_p):
        #   Formula (5.10) from NeST.pdf
        r"""
        Returns the reference pairing field for protons in uniform matter.
        This is an extrapolation between :math:`\Delta_{\mathrm{SM}}` and
        :math:`\Delta_{\mathrm{NeuM}}`. In limits :math:`\eta \rightarrow 0` reproduces
        symmetric matter and :math:`\eta \rightarrow 1`, the neutron matter.

        .. math::

            \Delta_p(\rho_n,\rho_p) =
            \Delta_{\mathrm{SM}}(\rho_n+\rho_p)
            \left( 1 - |\eta| \right)
            - \Delta_{\mathrm{NeuM}}(\rho_n) \eta \frac{\rho_p}{\rho_n+\rho_p}


        Args:
            rho_n (float): neutron density :math:`\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: pairing field for protons :math:`\Delta_p` [MeV]

        See also:
            :func:`.neutron_pairing_field`
            :func:`.symmetric_pairing_field`
            :func:`.neutron_ref_pairing_field`
        """
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        rho, eta = rhoEta(rho_n, rho_p) #eta = rho_n - rho_p
        rho = np.asarray(rho + DENSEPSILON, dtype=float)
        return _clamp_nonnegative(self.symmetric_pairing_field(rho_n, rho_p)*(1-np.abs(eta/rho))
            -self.neutron_pairing_field(rho_p)*rho_p/rho*eta/rho)


    # ================================
    #   Alternative interpolation schemes (delta_n, delta_p from delta_sm, delta_neu)
    # ================================
    # TODO(Adarsh): neutron_ref_pairing_field() / proton_ref_pairing_field() above
    # implement one interpolation scheme (Formula 5.10 from NeST.pdf) for getting
    # delta_n and delta_p from delta_sm (symmetric_pairing_field) and delta_neu
    # (neutron_pairing_field). Add two more schemes here, following:
    # https://link.springer.com/article/10.1140/epja/s10050-025-01503-x

    def neutron_ref_pairing_field_eq2(self, rho_n, rho_p):
        r"""
        Returns the reference pairing field for neutrons in uniform matter
        using the interpolation scheme from Eq. (2) of
        https://link.springer.com/article/10.1140/epja/s10050-025-01503-x

        This interpolates between symmetric matter and neutron matter,
        with asymmetry dependence governed by :math:`\delta`.

        Used in all Bsk parameterizations since Bsk17. Currently in use in WBSK.

        .. math::

            \Delta_n(\rho_n,\rho_p) =
            \Delta_{\mathrm{SM}}(\rho)
            \left(1 - |\delta|\right)
            + \delta \frac{\rho_n}{\rho}
            \Delta_{\mathrm{NeuM}}(\rho_n)

        where :math:`\rho = \rho_n + \rho_p` and
        :math:`\delta = (\rho_n - \rho_p)/\rho`.

        Args:
            rho_n (float): neutron density :math:`\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: pairing field for neutrons :math:`\Delta_n` [MeV]

        See also:
            :func:`.proton_ref_pairing_field_eq2`
            :func:`.symmetric_pairing_field`
            :func:`.neutron_pairing_field`
        """

        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        rho, eta = rhoEta(rho_n, rho_p) #eta = rho_n - rho_p
        rho = np.asarray(rho + DENSEPSILON, dtype=float)
        delta = eta/rho
        pairing = ((1-np.abs(delta))*self.symmetric_pairing_field(rho_n,rho_p)) + (delta* (rho_n/rho) * self.neutron_pairing_field(rho_n))
        return _clamp_nonnegative(pairing)

    def proton_ref_pairing_field_eq2(self, rho_n, rho_p):

        r"""
        Returns the reference pairing field for protons in uniform matter
        using the interpolation scheme from Eq. (2) of
        https://link.springer.com/article/10.1140/epja/s10050-025-01503-x

        This interpolates between symmetric matter and neutron matter,
        with asymmetry dependence governed by :math:`\delta`.

        Used in all Bsk parameterizations since Bsk17. Currently in use in WBSK.

        .. math::

            \Delta_p(\rho_n,\rho_p) =
            \Delta_{\mathrm{SM}}(\rho)
            \left(1 - |\delta|\right)
            - \delta \frac{\rho_p}{\rho}
            \Delta_{\mathrm{NeuM}}(\rho_p)

        where :math:`\rho = \rho_n + \rho_p` and
        :math:`\delta = (\rho_n - \rho_p)/\rho`.

        Args:
            rho_n (float): neutron density :math:`\rho_n` [fm :sup:`-3`]
            rho_p (float): proton density :math:`\rho_p` [fm :sup:`-3`]

        Returns:
            float: pairing field for protons :math:`\Delta_p` [MeV]

        See also:
            :func:`.neutron_ref_pairing_field_eq2`
            :func:`.symmetric_pairing_field`
            :func:`.neutron_pairing_field`
        """
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        rho, eta = rhoEta(rho_n, rho_p) #eta = rho_n - rho_p
        rho = np.asarray(rho + DENSEPSILON, dtype=float)
        delta = eta/rho
        pairing = ((1-np.abs(delta))*self.symmetric_pairing_field(rho_n,rho_p)) - (delta* (rho_p/rho) * self.neutron_pairing_field(rho_p))
        return _clamp_nonnegative(pairing)

    def neutron_ref_pairing_field_eq3(self, rho_n, rho_p):

        r"""
        Returns the reference pairing field for neutrons in uniform matter
        using the interpolation scheme from Eq. (3) of
        https://link.springer.com/article/10.1140/epja/s10050-025-01503-x

        This interpolates between symmetric matter and neutron matter,
        with asymmetry dependence governed by :math:`\delta`.

        Used in BskG3. 

        .. math::

            \Delta_n(\rho_n,\rho_p) =
            \Delta_{\mathrm{SM}}(\rho)
            \left(1 - |\delta|\right)
            + |\delta| \,
            \Delta_{\mathrm{NeuM}}(\rho_n)

        where :math:`\delta = (\rho_n - \rho_p)/\rho`.

        Args:
            rho_n (float): neutron density :math:`\rho_n` [fm :sup:`-3`]
            rho_p (float): proton density :math:`\rho_p` [fm :sup:`-3`]

        Returns:
            float: pairing field for neutrons :math:`\Delta_n` [MeV]

        See also:
            :func:`.proton_ref_pairing_field_eq3`
            :func:`.symmetric_pairing_field`
            :func:`.neutron_pairing_field`
        """
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        rho, eta = rhoEta(rho_n, rho_p) #eta = rho_n - rho_p
        rho = np.asarray(rho + DENSEPSILON, dtype=float)
        delta = eta/rho
        pairing = ((1-np.abs(delta))*self.symmetric_pairing_field(rho_n,rho_p)) + np.abs(delta) * self.neutron_pairing_field(rho_n)
        return _clamp_nonnegative(pairing)

    def proton_ref_pairing_field_eq3(self, rho_n, rho_p):

        r"""
        Returns the reference pairing field for protons in uniform matter
        using the interpolation scheme from Eq. (3) of
        https://link.springer.com/article/10.1140/epja/s10050-025-01503-x

        This interpolates between symmetric matter and neutron matter,
        with asymmetry dependence governed by :math:`\delta`.

        Used in BskG3.

        .. math::

            \Delta_p(\rho_n,\rho_p) =
            \Delta_{\mathrm{SM}}(\rho)
            \left(1 - |\delta|\right)
            + |\delta| \,
            \Delta_{\mathrm{NeuM}}(\rho_p)

        where :math:`\delta = (\rho_n - \rho_p)/\rho`.

        Args:
            rho_n (float): neutron density :math:`\rho_n` [fm :sup:`-3`]
            rho_p (float): proton density :math:`\rho_p` [fm :sup:`-3`]

        Returns:
            float: pairing field for protons :math:`\Delta_p` [MeV]

        See also:
            :func:`.neutron_ref_pairing_field_eq3`
            :func:`.symmetric_pairing_field`
            :func:`.neutron_pairing_field`
        """
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        rho, eta = rhoEta(rho_n, rho_p) #eta = rho_n - rho_p
        rho = np.asarray(rho + DENSEPSILON, dtype=float)
        delta = eta/rho
        pairing = ((1-np.abs(delta))*self.symmetric_pairing_field(rho_n,rho_p)) + np.abs(delta) * self.neutron_pairing_field(rho_p)
        return _clamp_nonnegative(pairing)

    def neutron_ref_pairing_field_eq6(self, rho_n, rho_p):

        r"""
        Returns the reference pairing field for neutrons in uniform matter
        using the interpolation scheme from Eq. (6) of
        https://link.springer.com/article/10.1140/epja/s10050-025-01503-x

        This interpolates between symmetric matter and neutron matter,
        with asymmetry dependence governed by :math:`\delta`.

        Used in BskG4.

        .. math::

            \Delta_n(\rho_n,\rho_p) =
            \Delta_{\mathrm{NeuM}}(\rho_n)
            \left[
            \frac{\Delta_{\mathrm{SM}}(\rho)}
                 {\Delta_{\mathrm{NeuM}}(\rho/2)}
            \right]^{(1 - \delta)}

        where :math:`\delta = (\rho_n - \rho_p)/\rho`.

        Args:
            rho_n (float): neutron density :math:`\rho_n` [fm :sup:`-3`]
            rho_p (float): proton density :math:`\rho_p` [fm :sup:`-3`]

        Returns:
            float: pairing field for neutrons :math:`\Delta_n` [MeV]

        See also:
            :func:`.proton_ref_pairing_field_eq6`
            :func:`.symmetric_pairing_field`
            :func:`.neutron_pairing_field`
        """
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        rho, eta = rhoEta(rho_n, rho_p) #eta = rho_n - rho_p
        rho = np.asarray(rho + DENSEPSILON, dtype=float)
        delta = eta/rho
        pairing = self.neutron_pairing_field(rho_n)* ((self.symmetric_pairing_field(rho_n,rho_p)/self.neutron_pairing_field(rho/2))**(1-delta))
        return _clamp_nonnegative(pairing)

    def proton_ref_pairing_field_eq6(self, rho_n, rho_p):

        r"""
        Returns the reference pairing field for protons in uniform matter
        using the interpolation scheme from Eq. (6) of
        https://link.springer.com/article/10.1140/epja/s10050-025-01503-x

        This interpolates between symmetric matter and neutron matter,
        with asymmetry dependence governed by :math:`\delta`.

        Used in BskG4.

        .. math::

            \Delta_p(\rho_n,\rho_p) =
            \Delta_{\mathrm{NeuM}}(\rho_p)
            \left[
            \frac{\Delta_{\mathrm{SM}}(\rho)}
                 {\Delta_{\mathrm{NeuM}}(\rho/2)}
            \right]^{(1 + \delta)}

        where :math:`\delta = (\rho_n - \rho_p)/\rho`.

        Args:
            rho_n (float): neutron density :math:`\rho_n` [fm :sup:`-3`]
            rho_p (float): proton density :math:`\rho_p` [fm :sup:`-3`]

        Returns:
            float: pairing field for protons :math:`\Delta_p` [MeV]

        See also:
            :func:`.neutron_ref_pairing_field_eq6`
            :func:`.symmetric_pairing_field`
            :func:`.neutron_pairing_field`
        """
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        rho, eta = rhoEta(rho_n, rho_p) #eta = rho_n - rho_p
        rho = np.asarray(rho + DENSEPSILON, dtype=float)
        delta = eta/rho
        pairing = self.neutron_pairing_field(rho_p)* ((self.symmetric_pairing_field(rho_n,rho_p)/self.neutron_pairing_field(rho/2))**(1+delta))
        return _clamp_nonnegative(pairing)


    # ================================
    #        Thermodynamics
    # ================================
    # Formulas from https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804:
    # Code formulas: A13, and dependent (A14, A15)
    def energy_per_nucleon(self, rho_n, rho_p):
        """
        Returns the energy per nucleon on infinite nuclear matter of given
        density of protons and neutrons, rho_p and rho_n, respectively, in MeV.
        Formula (A13) from https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804
        :cite:p:`chamel2009further`.

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: energy per nucleon :math:`e_n` [MeV]
        """
        p = self.params
        rho = rho_n+rho_p+DENSEPSILON
        kF = rho2kf(0.5*rho) # Formula (A14) from the paper is using rho/2
        eta = (rho_n-rho_p)/rho
        F_x_5 = 0.5*((1+eta)**(5./3.)+(1-eta)**(5./3.)) # Formula (A15)
        F_x_8 = 0.5*((1+eta)**(8./3.)+(1-eta)**(8./3.)) # Formula (A15)

        return (3*(HBARC**2)/20*(kF**2)*((np.power((1+eta),5/3)/MN
                                                    +np.power((1-eta),5/3)/MP))
                + p.t0/8.*rho*(3-(2*p.x0 + 1)*eta**2)
                + 3.*p.t1/40*rho*(kF**2)*((2+p.x1)*F_x_5 -(0.5+p.x1)*F_x_8)
                + 3./40*((2.*p.t2+p.t2x2)*F_x_5 +(0.5*p.t2+p.t2x2)*F_x_8)*rho*(kF**2)
                + p.t3/48*np.power(rho,(p.alpha+1))*(3-(1+2*p.x3)*(eta**2))
                + 3.*p.t4/40.*(kF**2)*np.power(rho,p.beta+1)*((2+p.x4)*F_x_5-(0.5+p.x4)*F_x_8)
                + 3.*p.t5/40.*(kF**2)*np.power(rho,p.gamma+1)*((2+p.x5)*F_x_5+(0.5+p.x5)*F_x_8))

    def testMe(self, rho):
        p = self.params
        rho = rho+DENSEPSILON
        kF = rho2kf(0.5*rho) # Formula (A14) from the paper is using rho/2
        r = 3*p.hbar2m_n/5*kF**2
        r =r+ 3./8*p.t0*rho
        r =r+ 3./80*(3*p.t1 + 5*p.t2 + 4*p.t2x2) * rho * kF**2
        r =r+ 1./16*p.t3*np.power(rho, p.alpha+1)
        r =r+ 9./80*p.t4*np.power(rho,p.beta+1)*kF**2
        r =r+ 3./80*p.t5*(5+4*p.x5)*np.power(rho,p.gamma+1)*kF**2
        return r

    def numerical_derivative_energy_per_nucleon_n(self, rho_n):
        """
        First derivative of energy per nucleon, calculated numerically using
        :func:`numpy.gradient`.

        Args:
        rho_n (float): rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns
            float: derivative of energy per neutron, :math:`\\frac{\\delta e_n}{\\delta \\rho}` [MeV]

        See also:
            :func:`.energy_per_nucleon`

        """
        e = self.energy_per_nucleon(rho_n, 0.)
        return np.gradient(e, rho_n)

    def numerical_second_derivative_energy_per_nucleon_n(self, rho_n):
        """
        Second derivative of energy per nucleon, calculated numerically using
        :func:`numpy.gradient`.

        Args:
        rho_n (float): rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns
            float: second derivative of energy per neutron, :math:`\\frac{\\delta^2 e_n}{\\delta \\rho^2}` [MeV]

        See also:
            :func:`.energy_per_nucleon`

        """
        e = self.numerical_derivative_energy_per_nucleon_n(rho_n)
        return np.gradient(e, rho_n)

    def numerical_derivative_epsilon_n(self, rho_n):
        """
        Derivative of energy density :math:`\\epsilon` with respect to density, calculated numerically using
        :func:`numpy.gradient` for neutron matter and using the :func:`.energy_per_nucleon` function, setting proton
        density :math:`\\rho_p` to 0.

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: Derivative of :math:`\\epsilon`

        See also:
            :func:`.energy_per_nucleon`
        """
        epsilon_neum = rho_n*(self.energy_per_nucleon(rho_n, 0.)+MN)
        return np.gradient(epsilon_neum, rho_n)

    def numerical_pressure_n(self, rho_n):
        """
        Pressure :math:`P` in neutron matter, calculated using the :func:`numpy.gradient`
        function to calculate the derivative of data calculated by :func:`.energy_per_nucleon`.

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: Pressure :math:`P`

        See also:
            :func:`.energy_per_nucleon`
        """
        derivative = np.gradient(self.energy_per_nucleon(rho_n, 0.), rho_n)
        return rho_n**2*derivative

    def numerical_derivative_pressure_n(self, rho_n):
        """
        Derivative of pressure :math:`P` with respect to density for neutron matter,
        calculated numerically using :func:`numpy.gradient`.

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: Derivative of :math:`P`

        See also:
            :func:`.numerical_pressure_n`
        """
        return np.gradient(self.pressure_n(rho_n), rho_n)

    def numerical_speed_of_sound_n(self, rho_n):
        """
        Velocity of sound in neutron matter dependent on total matter density equal
        to neutron density :math:`\\rho_n`, given as a percentage of speed of light.
        Calculated based on numerical data.

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: velocity of sound :math:`v_s`

        See also:
            :func:`.numerical_derivative_pressure_n`
            :func:`.numerical_derivative_epsilon_n`
        """
        return np.sqrt(self.numerical_derivative_pressure_n(rho_n)/self.numerical_derivative_epsilon_n(rho_n))*100

    #ALGEBRAIC ATTEMPT
    def energy_per_nucleon_n(self, rho):
        """
        Energy per nucleon on infinite neutron matter of density :math:`\\rho_n`.
        Derived from formula (A13) from https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804
        :cite:p:`chamel2009further`.

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: energy per neutron :math:`e_n` [MeV]
        """
        p = self.params
        c = (3*np.pi**2)**(2./3.)
        return (3./5.*p.hbar2m_n*c*rho**(2./3.)
                +1./4.*p.t0*(1-p.x0)*rho
                +3./40.*p.t1*(1-p.x1)*c*rho**(5./3.)
                +9./40.*(p.t2+p.t2x2)*c*rho**(5./3.)
                +1./24.*p.t3*(1-p.x3)*rho**(p.alpha+1)
                +3./40.*p.t4*(1-p.x4)*c*rho**(p.beta+5./3.)
                +9./40.*p.t5*(1+p.x5)*c*rho**(p.gamma+5./3.)
                )

    def derivative_energy_per_nucleon_n(self, rho):
        """
        First derivative of energy per nucleon (for neutrons only).

        Args:
            rho_n (float): rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: derivative of energy per neutron, :math:`\\frac{\\delta e_n}{\\delta \\rho}` [MeV]

        See also:
            :func:`.energy_per_nucleon_n`
        """
        p = self.params
        c = (3*np.pi**2)**(2./3.)
        return (2./5.*p.hbar2m_n*c*rho**(-1./3.)
                +1./4.*p.t0*(1-p.x0)
                +1./8.*p.t1*(1-p.x1)*c*rho**(2./3.)
                +3./8.*(p.t2+p.t2x2)*c*rho**(2./3.)
                +1./24.*p.t3*(1-p.x3)*(p.alpha+1)*rho**(p.alpha)
                +3./40.*p.t4*(1-p.x4)*c*(p.beta+5./3.)*rho**(p.beta+2./3.)
                +9./40.*p.t5*(1+p.x5)*c*(p.gamma+5./3.)*rho**(p.gamma+2./3.)
                )

    def second_derivative_energy_per_nucleon_n(self, rho):
        """
        Second derivative of energy per nucleon (for neutrons only).

        Args:
            rho_n (float): rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: second derivative of energy per neutron, :math:`\\frac{\\delta^2 e_n}{\\delta \\rho^2}` [MeV]

        See also:
            :func:`.energy_per_nucleon_n`
            :func:`.derivative_energy_per_nucleon_n`
        """
        p = self.params
        c = (3*np.pi**2)**(2./3.)
        return (-2./15.*p.hbar2m_n*c*rho**(-4./3.)
                +1./12.*p.t1*(1-p.x1)*c*rho**(-1./3.)
                +1./4.*(p.t2+p.t2x2)*c*rho**(-1./3.)
                +1./24.*p.t3*(1-p.x3)*(p.alpha+1)*p.alpha*rho**(p.alpha-1)
                +3./40.*p.t4*(1-p.x4)*c*(p.beta+5./3.)*(p.beta+2./3.)*rho**(p.beta-1./3.)
                +9./40.*p.t5*(1+p.x5)*c*(p.gamma+5./3.)*(p.gamma+2./3.)*rho**(p.gamma-1./3.)
                )

    def pressure_n(self, rho):
        """
        Pressure :math:`P` in neutron matter.

        Args:
            rho (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: Pressure :math:`P`

        See also:
            :func:`.derivative_energy_per_nucleon_n`
        """
        return rho**2*self.derivative_energy_per_nucleon_n(rho)

    def pressure_derivative_n(self, rho):
        """
        Calculates the derivative of pressure :math:`\\frac{\\delta P}{\\delta \\rho}` in neutron matter.

        Args:
            rho (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: derivative of pressure :math:`\\frac{\\delta P}{\\delta \\rho}`

        See also:
            :func:`.derivative_energy_per_nucleon_n`
            :func:`.second_derivative_energy_per_nucleon_n`
        """
        return 2*rho*self.derivative_energy_per_nucleon_n(rho) + rho**2*self.second_derivative_energy_per_nucleon_n(rho)

    def epsilon_derivative_n(self, rho):
        """
        Calculates the derivative of energy density :math:`\\frac{\\delta \\epsilon}{\\delta \\rho}` in neutron matter.

        Args:
            rho (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: derivative of energy density :math:`\\frac{\\delta \\epsilon}{\\delta \\rho}`

        See also:
            :func:`.derivative_energy_per_nucleon_n`
            :func:`.energy_per_nucleon`
        """
        return MN + self.energy_per_nucleon(rho, 0.) + rho*self.derivative_energy_per_nucleon_n(rho)

    def speed_of_sound_n(self, rho_n):
        """
        Velocity of sound in neutron matter, dependent on total matter density equal
        to neutron density :math:`\\rho_n`, given as a percentage of speed of light.

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: velocity of sound :math:`v_s`

        See also:
            :func:`.pressure_derivative_n`
            :func:`.epsilon_derivative_n`
        """
        return np.sqrt(self.pressure_derivative_n(rho_n)/self.epsilon_derivative_n(rho_n))*100


    # ================================
    #       Effective masses
    # ================================
    # Implementation based on:
    # https://journals.aps.org/prc/pdf/10.1103/PhysRevC.82.035804
    # Equation 10, plotted as function of rho (see Fig 5, inset)
    # 
    # NOTE: Functions calculate effective masses for both species (q=n or p)
    # where total density rho = rho_n + rho_p
    #
    # Common use cases:
    # - Neutron matter (NeuM): rho = rho_n, rho_p = 0
    # - Symmetric matter (SM): rho = 2*rho_n, rho_p = rho_n
    #
    # TODO: Add plotting functions to visualize M* vs rho for different asymmetries

    def isoscalarM(self, rho_n, rho_p):
        """
        Calculates effective isoscalar mass M_s for a given uniform system
        with neutron and proton densities rho_n and rho_p respectively.

        .. math::

            \\frac{M^*_s}{M} = 2 \\left( \\frac{M}{M_n^*} + \\frac{M}{M_p^*} \\right)^{-1}

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: effective isoscalar mass in units of the bare mass, :math:`M_{s}^{*}/M` [1]

        See also:
            :func:`.effMn`
            :func:`.effMp`

        """
        return 2/((1./self.effMn(rho_n, rho_p))+(1./self.effMp(rho_n, rho_p)))

    def isovectorM(self, rho_n, rho_p):
        """
        Calculates effective isovector mass :math:`M_v` for a given uniform system
        with neutron and proton densities rho_n, rho_p respectively.

        .. math::

            M^{*}_v = M^{*}_s M_p^{*} \\frac{ 2 \\rho_p-\\rho}{2 M_p^*\\rho_p-\\rho}

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: effective isovector mass in units of the bare mass, :math:`M_{v}^{*}/M` [1]

        See also:
            :func:`.effMn`
            :func:`.effMp`
            :func:`.isoscalarM`
        """
        Mp = self.effMp(rho_n, rho_p)
        rho = rho_n + rho_p + DENSEPSILON
        return self.isoscalarM(rho_n, rho_p)*Mp*(2*rho_p-rho)/(2*Mp*rho_p-rho)

    def effMn(self, rho_n, rho_p):
        """
        Returns the effective mass of a neutron in nuclear medium.

        .. math::

            \\frac{M_n^*}{M_n} = \\frac{\\hbar^2}{2 M_n} \\frac{1}{B_n(\\rho_n,\\rho_p)}

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: neutron effective mass in units of the bare neutron mass, :math:`M_{n}^{*}/M_n` [1]

        See also:
            :func:`.B_q`
        """
        p = self.params
        return p.hbar2m_n/self.B_q(rho_n, rho_p,'n')

    def effMp(self, rho_n, rho_p):
        """
        Returns the effective mass of a proton in nuclear medium.

        .. math::

            \\frac{M_p^*}{M_p} = \\frac{\\hbar^2}{2 M_p} \\frac{1}{B_p(\\rho_n,\\rho_p)}

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: proton effective mass in units of the bare proton mass, :math:`M_{p}^{*}/M_p` [1]

        See also:
            :func:`.B_q`
        """
        p = self.params
        return p.hbar2m_p/self.B_q(rho_n, rho_p,'p')

    def mu_q(self, rho_n, rho_p, q):
        # Eq. taken from S. Goriely, N. Chamel, and J. M. Pearson, Phys. Rev. Lett. 102, 152503 (2009)
        """
        Calculates the chemical potential :math:`\\mu_q` defined with the wavevector
        :math:`k_{F,q}` and the effective mass :math:`M^*_q` :cite:`chamel2009pairing`.

        .. math::

            \\mu_q = \\frac{\\hbar^2 k_{F,q}^2}{2M^*_q} = B_q k_{F,q}^2

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
            q (str): nucleon type choice ('p' - proton, or 'n' - neutron)

        Returns:
             float: chemical potential :math:`\\mu` [MeV]
        """

        rho_n = np.asarray(rho_n)
        rho_p = np.asarray(rho_p)
        if(q=='n'):
            rho = rho_n
        elif(q=='p'):
            rho = rho_p
        else:
            sys.exit('# ERROR: Nucleon q must be either n or p')
        mu_q = self.B_q(rho_n, rho_p, q)*rho2kf(rho)**2
        mu_q = np.asarray(mu_q)
        # Handle assignment for both scalar and array
        if mu_q.shape == ():
            if mu_q == 0:
                mu_q = NUMZERO
            if rho == 0:
                mu_q = NUMZERO
            return float(mu_q)
        else:
            # boolean masks, broadcast: rho may be a scalar when the other density is an array
            mu_q[mu_q == 0] = NUMZERO
            mu_q[np.broadcast_to(rho == 0, mu_q.shape)] = NUMZERO
            return mu_q

    def U_q(self, rho_n, rho_p,q):
        #     Formula (5.14) from NeST.pdf
        """
        Returns the mean field potential from density :math:`\\rho` variation.

        rho_q is either rho_n or rho_p.

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
            q (str): nucleon type choice ('p' - proton, or 'n' - neutron)

        Returns:
            float: Mean field potential :math:`U_q` [MeV]
        """
        p = self.params
        if(q=='n'):
            rho_q = rho_n + DENSEPSILON
            rho_q_prime = rho_p + DENSEPSILON
        elif(q=='p'):
            rho_q = rho_p + DENSEPSILON
            rho_q_prime = rho_n + DENSEPSILON
        else:
            sys.exit('# ERROR: Nucleon q must be either n or p')
        rho = rho_n + rho_p + DENSEPSILON
        return p.t0*((1+0.5*p.x0)*rho-(0.5+p.x0)*rho_q)+p.t3/12.*np.power(rho,(p.alpha-1))*(
            ((0.5+0.5*p.x3)*rho**2*(p.alpha+2))-(0.5+p.x3)*(2*rho*rho_q*p.alpha*
                                                      (rho_q_prime)**2))

    def B_q(self, rho_n, rho_p, q):
        #    Formula (5.13) from NeST.pdf
        """
        Returns the mean field potential B_q (coming from variation over kinetic
        density, or effective mass).

        .. math::

        	B_q  = \\frac{\\hbar^2}{2 M^*_q}
        	 = \\frac{\\hbar^2}{2 M_q}  + C^\\tau_0\\rho  + C^\\tau_1 (\\rho_{q} - \\rho_{q'})


        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
            q (str): nucleon type choice ('p' - proton, or 'n' - neutron)

        Returns:
            float: :math:`B_q = \\hbar^2/(2 M^*_q)` [MeV fm :sup:`2`]
        """
        p = self.params
        if(q=='n'):
            rho_q = rho_n + DENSEPSILON
            HBAR2M_q = p.hbar2m_n
        elif(q=='p'):
            rho_q = rho_p + DENSEPSILON
            HBAR2M_q = p.hbar2m_p
        else:
            sys.exit('# ERROR: Nucleon q must be either n or p')
        rho = rho_n + rho_p
        return (HBAR2M_q
                 + p.t1/4.*((1.+p.x1/2.)*rho - (1./2.+p.x1)*rho_q)
                 + p.t4/4.*((1.+p.x4/2.)*rho - (1./2.+p.x4)*rho_q)*np.power(rho,p.beta)
                 + 1./4.*((p.t2+p.t2x2/2.)*rho+(1./2.*p.t2+p.t2x2)*rho_q)
                 + p.t5/4.*((1.+p.x5/2.)*rho + (1./2.+p.x5)*rho_q)*np.power(rho,p.gamma)
                 )


    # ================================
    #        Energy functional
    # ================================
    # Implementation based on:
    # PHYSICAL REVIEW C 104, 055801 (2021)
    # Formulas: 9-14, 23, A8-A10
    # 
    # NOTE: Current implementation uses uniform matter approximation.
    # Future extensions will include:
    # - Non-uniform densities (TAU, MU, J) from simulation data
    # - Gradient corrections
    # - Spin-orbit and tensor terms
    #
    # TODO: Extend to non-uniform matter with kinetic densities from data files


    def C0_rho(self, rho_n, rho_p):
        """
        Calculates the energy functional :math:`C^{\\rho}_0` coefficient :cite:`bender2002gamow`.

        .. math::

        	 C^0_\\rho(\\rho) = \\frac{3}{8} t_0 + \\frac{3}{48} t_3 \\rho^\\alpha

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: :math:`C^{\\rho}_0` coefficient
        """
        p = self.params
        return 3. / 8. * p.t0 + 3. / 48. * p.t3 * np.power(rho_n + rho_p, p.alpha)


    def C1_rho(self, rho_n, rho_p):
        """
        Calculates the energy functional :math:`C^{\\rho}_1` coefficient :cite:`bender2002gamow`.

        .. math::

            C^1_\\rho(\\rho) = -\\frac{1}{4} t_0 \\left(\\frac{1}{2} + x_0\\right) - \\frac{1}{24} t_3 \\left(\\frac{1}{2} + x_3 \\right) \\rho^\\alpha

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: :math:`C^{\\rho}_1` coefficient
        """
        p = self.params
        return -1. / 4. * p.t0 * (1. / 2. + p.x0) - 1. / 24. * p.t3 * (1. / 2. + p.x3) * np.power(rho_n + rho_p, p.alpha)


    def C0_tau(self, rho_n, rho_p):
        """
        Calculates the energy functional :math:`C^{\\tau}_0` coefficient :cite:`bender2002gamow`.

        .. math::

        	 C^0_\\tau(\\rho) =  \\frac{3}{16} t_1 + \\frac{1}{4} t_2 \\left(\\frac{5}{4} + x_2\\right) + \\frac{3}{16} t_4 \\rho^\\beta + \\frac{1}{4} t_5 \\left(\\frac{5}{4} + x_5 \\right) \\rho^\\gamma

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: :math:`C^{\\tau}_0` coefficient
        """
        p = self.params
        rho = rho_n + rho_p
        return (3. / 16.) * p.t1 + 5. / 16. * p.t2 + 1. / 4. * p.t2x2 + 3. / 16. * p.t4 * np.power(rho, p.beta) + 1. / 4. * p.t5 * (1.25 + p.x5) * np.power(rho, p.gamma)


    def C1_tau(self, rho_n, rho_p):
        """
        Calculates the energy functional :math:`C^{\\tau}_1` coefficient :cite:`bender2002gamow`.

        .. math::

        	 C^1_\\tau(\\rho) = -\\frac{1}{8} t_1 \\left( \\frac{1}{2} + x_1 \\right) + \\frac{1}{8} t_2 \\left(\\frac{1}{2} + x_2\\right) - \\frac{1}{8} t_4 \\rho^\\beta \\left( \\frac{1}{2} + x_4 \\right)  + \\frac{1}{8} t_5 \\left(\\frac{1}{2} + x_5 \\right) \\rho^\\gamma

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: :math:`C^{\\tau}_1` coefficient
        """
        p = self.params
        rho = rho_n + rho_p
        return (-1. / 8.) * p.t1 * (1. / 2. + p.x1) + 1. / 16. * p.t2 + 1. / 8. * p.t2x2 - 1. / 8. * p.t4 * np.power(rho, p.beta) * (1. / 2. + p.x4) + 1. / 8. * p.t5 * np.power(rho, p.gamma) * (1. / 2. + p.x5)

    def epsilon_rho_np(self, rho_n, rho_p):
        """
        Energy functional :math:`\\epsilon_{\\rho}` for particle matter, related to
        the interaction of the nucleons with the background matter density.

        .. math::

        	\\varepsilon_\\rho(\\rho_n,\\rho_p)
        	= C^\\rho_0(\\rho_n+\\rho_p) \\left( \\rho_n + \\rho_p \\right)^2
        	 + C^\\rho_1(\\rho_n+\\rho_p) \\left( \\rho_n - \\rho_p \\right)^2

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components

        Returns:
            float: energy functional :math:`\\epsilon_{\\rho}`
        """
        rho_n = np.asarray(rho_n)
        rho_p = np.asarray(rho_p)
        rho = rho_n + rho_p + DENSEPSILON
        return self.C0_rho(rho_n, rho_p) * np.power(rho, 2) + self.C1_rho(rho_n, rho_p) * np.power(rho_n - rho_p, 2)


    def epsilon_tau_np(self, rho_n, rho_p, tau_n, tau_p, jsum2, jdiff2):
        """
        Energy functional :math:`\\epsilon_{\\tau}` for particle matter,
        related to the density-dependent effective mass. It gives rise to current couplings.

        .. math::

        	\\varepsilon_\\tau(\\rho_n,\\tau_n,{ j}_n,\\rho_p,\\tau_p,{ j}_p)

            = C^\\tau_0(\\rho_n+\\rho_p) \\left[ (\\rho_n + \\rho_p)(\\tau_n + \\tau_p) - ({ j}_n + { j}_p)^2 \\right]

            + C^\\tau_1(\\rho_n+\\rho_p) \\left[ (\\rho_n - \\rho_p)(\\tau_n - \\tau_p) - ({ j}_n - { j}_p)^2 \\right].

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
            tau_n (float): kinetic density :math:`\\tau` [fm :sup:`-5`]
            tau_n (float): kinetic density :math:`\\tau` [fm :sup:`-5`]
            jsum2 (float): sum of momentum density/current vectors :math:`j` [fm :sup:`-3`]
            jdiff2 (float) : difference of momentum density/current vectors :math:`j` [fm :sup:`-3`]

        Returns:
            float: energy functional :math:`\\epsilon_{\\tau}`
        """
        rho_n = np.asarray(rho_n)
        rho_p = np.asarray(rho_p)
        tau_n = np.asarray(tau_n)
        tau_p = np.asarray(tau_p)
        jsum2 = np.asarray(jsum2)
        jdiff2 = np.asarray(jdiff2)
        return self.C0_tau(rho_n, rho_p) * ((rho_n + rho_p) * (tau_n + tau_p) - jsum2) + self.C1_tau(rho_n, rho_p) * ((rho_n - rho_p) * (tau_n - tau_p) - jdiff2)


    def epsilon_delta_rho_np(self, rho_n, rho_p, rho_grad_n_square, rho_grad_p_square, rho_grad_square):
        """
        Energy functional :math:`\\epsilon_{\\Delta \\rho}` for particle matter,
        related to the interaction of the nucleons with the background fluctuations.

        .. math::
           \\begin{equation*}
           \\begin{split}
        	\\varepsilon_{\\Delta\\rho}&(\\rho_n,\\vec\\nabla\\rho_n,\\rho_p,\\vec\\nabla\\rho_p)	=
        	 +\\frac{3}{16} t_1 \\left[
        	  \\left(          1 + \\frac{1}{2}x_1\\right)         \\left( \\nabla(\\rho_n+\\rho_p)   \\right)^2
        	 -\\left(\\frac{1}{2} +           x_1 \\right)  \\sum_q \\left( \\nabla\\rho_q \\right)^2
        	 \\right] \\\\
        	  &- \\frac{1}{16}t_2 \\left[
        	     \\left(           1 + \\frac{1}{2} x_2\\right)          \\left( \\nabla(\\rho_n+\\rho_p)  \\right)^2
        	   + \\left( \\frac{1}{2} +             x_2\\right) \\sum_{q} \\left( \\nabla\\rho_q\\right)^2  \\right] \\\\
        	  &+\\frac{3}{16} t_4  (\\rho_n+\\rho_p)^\\beta \\left[
        	  \\left(1 + \\frac{1}{2}x_4\\right)           \\left( \\nabla(\\rho_n+\\rho_p)   \\right)^2
        	 - \\left(\\frac{1}{2} + x_4 \\right)  \\sum_q  \\left( \\nabla\\rho_q \\right)^2
        	 \\right] \\\\
        	 &+ \\frac{\\beta}{8}t_4 (\\rho_n+\\rho_p)^{\\beta-1} \\left[
        	 \\left(1 + \\frac{1}{2}x_4\\right) (\\rho_n+\\rho_p)\\left(\\nabla(\\rho_n+\\rho_p)\\right)^2
        	 - \\left(\\frac{1}{2} + x_4 \\right)
        	 \\nabla(\\rho_n+\\rho_p)\\cdot\\left(\\sum_q \\rho_q\\nabla\\rho_q\\right)
        	 \\right] \\\\
        	 &- \\frac{1}{16}t_5 (\\rho_n+\\rho_p)^\\gamma \\left[
        	   \\left(          1 + \\frac{1}{2} x_5\\right) \\left( \\nabla(\\rho_n+\\rho_p)\\right)^2
        	  +\\left(\\frac{1}{2} +             x_5\\right) \\sum_{q} \\left( \\nabla\\rho_q\\right)^2
        	 \\right].
           \\end{split}
           \\end{equation*}


        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
            rho_grad_n_square (float): neutron density gradient squared :math:`\\nabla \\rho` [fm :sup:`-8`]
            rho_grad_p_square (float): proton density gradient squared :math:`\\nabla \\rho` [fm :sup:`-8`]
            rho_grad_square (float): total density gradient squared :math:`\\nabla \\rho` [fm :sup:`-8`]

        Returns:
            float: energy functional :math:`\\epsilon_{\\Delta \\rho}`
        """
        p = self.params
        rho_n = np.asarray(rho_n)
        rho_p = np.asarray(rho_p)
        rho_grad_n_square = np.asarray(rho_grad_n_square)
        rho_grad_p_square = np.asarray(rho_grad_p_square)
        rho_grad_square = np.asarray(rho_grad_square)
        rho = rho_n+rho_p + DENSEPSILON
        grad_rho_n_rho_p = 0.5*(rho_grad_p_square-rho_grad_n_square-rho_grad_p_square)

        return (3./16.*p.t1*((1.+0.5*p.x1)*rho_grad_square-(0.5+p.x1)*(rho_grad_n_square+rho_grad_p_square))
                -1./16.*p.t2x2*(0.5*rho_grad_square+rho_grad_n_square+rho_grad_p_square)
                -1./16.*p.t2*(rho_grad_square+0.5*rho_grad_n_square+0.5*rho_grad_p_square)
                +3./16.*p.t4*np.power(rho,p.beta)*((1.+0.5*p.x4)*rho_grad_square-(0.5+p.x4)*(rho_grad_n_square+rho_grad_p_square))
                -1./16.*p.t5*np.power(rho,p.gamma)*((1.+0.5*p.x5)*rho_grad_square+(0.5+p.x5)*(rho_grad_n_square+rho_grad_p_square))
                +p.beta/8.*p.t4*np.power(rho,p.beta-1.)*((1.+0.5*p.x4)*rho*rho_grad_square-(0.5+p.x4)
                                          *(rho_n*rho_grad_n_square+rho_p*rho_grad_p_square+rho*grad_rho_n_rho_p))
                )


    def epsilon_np(self, rho_n, rho_p, rho_grad_n, rho_grad_p, tau_n, tau_p, jsum2, jdiff2, nu_n, nu_p, kappa_n, kappa_p):
        """
        Calculates the total energy functional :math:`\\epsilon`, with a limit for
        single-particle energies :math:`\\epsilon_{\\Lambda}` = 6.5 MeV.

        .. math::
            \\varepsilon(\\rho_n,\\vec\\nabla\\rho_n, \\tilde{\\rho}_n,\\tau_n,{ j}_n,\\rho_p,\\vec\\nabla\\rho_p,\\tilde{\\rho}_p,\\tau_p,{ j}_p) = \\nonumber \\\\
         = \\frac{\\hbar^2}{2 M_n} \\tau_n + \\frac{\\hbar^2}{2 M_p} \\tau_p
          + \\varepsilon_\\rho(\\rho_n,\\rho_p)
          + \\varepsilon_\\tau(\\rho_n,\\tau_n,{ j}_n,\\rho_p,\\tau_p,{ j}_p) \\nonumber \\\\
         + \\varepsilon_{\\Delta\\rho}(\\rho_n,\\vec\\nabla\\rho_n,\\rho_p,\\vec\\nabla\\rho_p)
         + \\varepsilon_\\pi(\\rho_n,\\vec\\nabla\\rho_n,\\tilde{\\rho}_n,\\rho_p,\\vec\\nabla\\rho_p,\\tilde{\\rho}_p),

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
            rho_grad_n (float): neutron density gradient :math:`\\nabla \\rho` [fm :sup:`-4`]
            rho_grad_p (float): proton density gradient :math:`\\nabla \\rho` [fm :sup:`-4`]
            tau_n (float): kinetic density :math:`\\tau` [fm :sup:`-5`]
            tau_n (float): kinetic density :math:`\\tau` [fm :sup:`-5`]
            jsum2 (float): sum of momentum density/current vectors :math:`j` [fm :sup:`-3`]
            jdiff2 (float) : difference of momentum density/current vectors :math:`j` [fm :sup:`-3`]
            nu_n (float): neutron anomalous density :math:`\\nu` [fm :sup:`-3`]
            nu_p (float): proton anomalous density :math:`\\nu` [fm :sup:`-3`]
            kappa_n (float):
            kappa_p (float):
                what is kappa? (no Eq.9 in Ref.41)

        TO DO: write eq for jsum2 and jdiff2

        Returns
            float: energy functional :math:`\\epsilon`
        """
        rho_n = np.asarray(rho_n)
        rho_p = np.asarray(rho_p)
        rho_grad_n = np.asarray(rho_grad_n)
        rho_grad_p = np.asarray(rho_grad_p)
        tau_n = np.asarray(tau_n)
        tau_p = np.asarray(tau_p)
        jsum2 = np.asarray(jsum2)
        jdiff2 = np.asarray(jdiff2)
        nu_n = np.asarray(nu_n)
        nu_p = np.asarray(nu_p)
        kappa_n = np.asarray(kappa_n)
        kappa_p = np.asarray(kappa_p)
        rho_grad_n_square = (rho_grad_n)**2
        rho_grad_p_square = (rho_grad_p)**2
        rho_grad_square = rho_grad_n_square + rho_grad_p_square
        return (HBARC**2/2./MN*tau_n + HBARC**2/2./MP*tau_p
                + self.epsilon_rho_np(rho_n, rho_p)
                + self.epsilon_delta_rho_np(rho_n, rho_p, rho_grad_n_square, rho_grad_p_square, rho_grad_square)
                + self.epsilon_tau_np(rho_n, rho_p, tau_n, tau_p, jsum2, jdiff2)
                + self.epsilon_pi_np(rho_n, rho_p, rho_grad_n, rho_grad_p, nu_n, nu_p, kappa_n, kappa_p)
                )


    def epsilon_pi_np(self, rho_n, rho_p, rho_grad_n, rho_grad_p, nu_n, nu_p, kappa_n,
                        kappa_p):
        """
        Energy functional :math:`\\epsilon_{\\pi}` for particle matter,
        related to pairing energy density.

        .. math::

        	\\varepsilon_\\pi(\\rho_n,\\vec\\nabla\\rho_n,\\tilde{\\rho}_n,\\rho_p,\\vec\\nabla\\rho_p,\\tilde{\\rho}_p)

            =\\frac{1}{4} f^\\pm_n \\left( v^{\\pi n}(\\rho_n,\\rho_p) + \\kappa_n|\\nabla\\rho_n|^2 \\right) \\tilde{\\rho_n}^2

            +\\frac{1}{4} f^\\pm_p \\left( v^{\\pi p}(\\rho_n,\\rho_p) + \\kappa_p|\\nabla\\rho_p|^2 \\right) \\tilde{\\rho_p}^2,

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
            rho_grad_n (float): neutron density gradient :math:`\\nabla \\rho_n` [fm :sup:`-4`]
            rho_grad_p (float): proton  density gradient :math:`\\nabla \\rho_p` [fm :sup:`-4`]
            nu_n (float): neutron anomalous density :math:`\\nu_n` [fm :sup:`-3`]
            nu_p (float): proton anomalous density :math:`\\nu_p` [fm :sup:`-3`]
            kappa_n: :math:`\\kappa_n` (float)
            kappa_p: :math:`\\kappa_p` (float)

        Returns:
            float: energy functional :math:`\\epsilon_{\\pi}` [MeV]

        See also:
            :func:`.v_pi`
        """
        return (1./4.*(self.v_pi(rho_n, rho_p, 'n') + kappa_n*(np.abs(rho_grad_n))**2)*nu_n**2
                +1./4.*(self.v_pi(rho_n, rho_p, 'p') + kappa_p*(np.abs(rho_grad_p))**2 )* nu_p**2)

    def v_pi(self, rho_n, rho_p, q):
        """
        Calculates pairing strength :math:`\\upsilon^{\\pi}_q` for neutrons
        or protons, for energies below 6.5 MeV.

        Based on Equation 14 from Phys Rev C 104.
        Check :cite:`pecak2021properties` and earlier papers of Chamel.

        .. note::
            This function handles both neutron (q='n') and proton (q='p') cases
            with the same implementation. Consider splitting into separate functions
            for better code organization if species-specific modifications are needed.

        .. math::

            v^{\\pi q}(\\rho_n,\\rho_p)
            = - \\frac{8 \\pi^2}{I_q(\\rho_n,\\rho_p) } B_q(\\rho_n,\\rho_p)^{3/2}

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
            q (str): nucleon type choice ('p' - proton, or 'n' - neutron)

        Returns:
            float: pairing strength :math:`\\upsilon^{pi}` [fm :sup:`-3`]

        See also:
            :func:`.I`
            :func:`.effMn`
            :func:`.effMp`
            :func:`.Lambda`
            :func:`.B_q`
        """
        #rho = rho_n + rho_p
        if(q=='n'):
            M = self.effMn(rho_n, rho_p)
        elif(q=='p'):
            M = self.effMp(rho_n, rho_p)
        else:
            sys.exit('# ERROR: Nucleon q must be either n or p')
        return -8.*np.pi**2/self.I(rho_n, rho_p, q) * (HBARC**2/2./M)**(3./2.)

    def I(self, rho_n, rho_p, q):
        # Equation 15 from Phys Rev C 104
        """
        Approximates the solution to the integral present in the original
        :math:`\\nu^{\\pi}` formula for single-particle energies below 6.5 MeV.
        Check :cite:`pecak2021properties` and earlier papers of Chamel.

        .. math::

            I_q(\\rho_n,\\rho_p)
            = \\sqrt{\\mu_q(\\rho_n,\\rho_p) } \\left( 2 \\ln{\\frac{2 \\mu_q(\\rho_n,\\rho_p) }{\\Delta_q(\\rho_n,\\rho_p) }} + \\Lambda\\left( \\frac{\\varepsilon_\\Lambda}{\\mu_q(\\rho_n,\\rho_p) } \\right) \\right)

        Args:
            rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components
            rho_p (float): proton density :math:`\\rho_p` [fm :sup:`-3`]; sum of both spin components
            q (str): nucleon type choice ('p' - proton, or 'n' - neutron)

        Returns:
            float: analytic integral solution

        See also:
            :func:`.Lambda`
            :func:`.mu_q`
            :func:`.neutron_ref_pairing_field`
            :func:`.proton_ref_pairing_field`
        """
        if(q=='n'):
            delta = self.neutron_ref_pairing_field(rho_n, rho_p)
        elif(q=='p'):
            delta = self.proton_ref_pairing_field(rho_n, rho_p)
        else:
            sys.exit('# ERROR: Nucleon q must be either n or p')

        mu = self.mu_q(rho_n, rho_p, q)
        I = np.asarray(mu)
        delta = np.asarray(delta)
        # Handle assignment for both scalar and array
        if I.shape == ():
            if delta == NUMZERO:
                I = NUMZERO
            elif delta > NUMZERO:
                I = np.sqrt(mu)*(2*np.log(2*mu/np.abs(delta))+Lambda(self.pairing.cutoff/mu))
            return float(I)
        else:
            x = np.where(delta==NUMZERO)
            I[x] = NUMZERO
            I = np.ma.masked_where(delta == 0, I)
            y = np.where(delta>NUMZERO)
            I[y] = np.sqrt(mu[y])*(2*np.log(2*mu[y]/np.abs(delta[y]))+Lambda(self.pairing.cutoff/mu[y]))
            return I
