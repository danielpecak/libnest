#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Info: who, where, when
# @author Daniel Pęcak <Daniel.Pecak@pw.edu.pl>
# Warsaw Technical University, Université Libre de Bruxelles
# On leave: Institute of Physics, Polish Academy of Sciences, Warsaw
# March 2022, Brussels
"""
Module: Units
=============
This module contains physical constants relevant for nuclear matter and neutron star physics. It also provides functions for converting one type of units into another.

Numerical
---------
Constants for dealing with extreme numerics:
 * ``NUMZERO`` = 1e-12 (the lowest numerical value allowed)
 * ``NUMINF``  = 1e30  (the largest numerical value allowed)
 * ``DENSEPSILON`` = 1e-12 this value is used to renormalize division and avoid dividing but extremely small numbers that cause numerical instabilities

.. note::
    ``DENSEPSILON`` is used whenever division by density has to be applied.
    Those cases usually will give zero in proper analytical limit.

Physical
--------
.. note::
    Physical constants are from NIST website:
     * https://physics.nist.gov/cgi-bin/cuu/Value?hbcmevf
     * https://physics.nist.gov/cgi-bin/cuu/Value?mnc2mev
     * https://physics.nist.gov/cgi-bin/cuu/Value?mpc2mev


========  =========== =================== ===========
Constant  Value       Unit                Description
========  =========== =================== ===========
HBARC     197.3269804 [MeV fm]            :math:`\hbar c`
MN        939.5654205 [MeV]               neutron mass :math:`m_n`
MP        938.2720882 [MeV]               proton  mass :math:`m_p`
HBAR2M_n  20.72124837 [MeV*fm :math:`^2`] :math:`\\hbar^2/(2 m_n)`
HBAR2M_p  20.74981092 [MeV*fm :math:`^2`] :math:`\\hbar^2/(2 m_p)`
========  =========== =================== ===========

Nucleon masses
--------------
The masses of proton, neutron and nucleon (which is their averaged mass). The nucleon mass is an average of proton and neutron mass:

.. math::
    m_N=\\frac{1}{2} (m_n + m_p) \\approx 1.67377585 \\cdot 10^{-27}.


..  csv-table::
    :header: "Mass", "Value"
    :widths: 10, 15

    neutron :math:`m_n`,    1.6749286 :math:`\cdot 10^{27}` kg
    proton  :math:`m_p`,    1.6726231 :math:`\cdot 10^{27}` kg
    nucleon :math:`m_N`,    1.6737759 :math:`\cdot 10^{27}` kg


Neutron star
------------
Some constants are relevant from the point of view of neutron star's physics.
One of them is the *neutron drip* density :math:`\\rho_{\\mathrm{ND}}` at which excesive neutrons are not bound to the nuclei anymore and form superfluid sea. This is how the border of outer crust and inner crust are defined.  Then, there is saturation density :math:`\\rho_0` at which the crust-core transition should occure. This is the density of nuclei.

..  csv-table::
    :header: "Variable", "Density", "[g cm :sup:`-3`]", "[fm :sup:`-3`]"
    :widths: 15, 15, 15, 15

    "RHOSAT", ":math:`\\rho_0`", "3 :math:`\\cdot 10^{14}`", 0.18
    "RHOND", ":math:`\\rho_{\\mathrm{ND}}`", "4 :math:`\\cdot 10^{11}`", 0.00042

List of functions
-----------------
"""
import numpy as np

DENSEPSILON = 1e-12
NUMZERO = 1e-12
NUMINF  = 1e30

# EXACT:
E  = 1.602176634e-19 # elementary charge [C]
kB = 1.380649e-23 # Boltzmann constant [J/K]
ALPHA=7.2973525693e-3 # fine-structure constants
HBARC=197.3269804  # \hbar c [MeV fm]
hbar22M0   =20.72  # neutron bare mass
MN   =939.56542052 # neutron mass [MeV]
MP   =938.27208816 # proton  mass [MeV]
HBAR2M_n = 20.721248369006936 # [MeV*fm<sup>2</sup>] 0.5*hbar^2/neutron mass
HBAR2M_p = 20.749810921501915 # [MeV*fm<sup>2</sup>] 0.5*hbar^2/proton mass


def KtoMev(temp):
    """
    Converts temperature :math:`T` units from Kelvins to energy :math:`E` units in MeVs by setting the Boltzmann constant to 1.

    .. math::
        E = \\frac{k_B}{1eV} 10^{-6} T \\approx 11604525006.1598 \\cdot T

    Args:
        temp (float): temperature :math:`T` [K]

    Returns:
        float: energy :math:`E` [MeV]

    See also:
        :func:`MeVtoK`
    """
    temp = np.asarray(temp)
    return temp/11604525006.1598

def MeVtoK(energy):
    """
    Converts energy :math:`E` units in MeVs to temperature :math:`T` units from Kelvins by setting the Boltzmann constant to 1.

    .. math::
        T = \\frac{e V}{k_B} 10^6 E \\approx 11604525006.1598 \\cdot E


    Args:
        energy (float): energy :math:`E` [MeV]

    Returns:
        float: temperature :math:`T` [K]

    See also:
        :func:`KtoMev`
    """
    energy = np.asarray(energy)
    return energy*11604525006.1598

def fm3togcm3(rho):
    """
    Function converts desnity units: fm :sup:`-3` into g cm :sup:`-3`. See more `here <#nucleon-masses>`_. The numerical factor

    .. math::
        \\frac{m_N}{\\mathrm{fm}^3} = 1.67377585 \\cdot  10^{15}\\frac{\\mathrm{g}}{\\mathrm{cm}^3}

    Args:
        rho (float): density :math:`\\rho` [fm :sup:`-3`]

    Returns:
        float: rho :math:`\\rho` [g cm :sup:`-3`]

    See also:
        :func:`gcm3tofm3`
    """
    rho = np.asarray(rho)
    return rho*1.67377585e15

def gcm3tofm3(rho):
    """
    Function converts desnity units: fm :sup:`-3` into g cm :sup:`-3`. See more `here <#nucleon-masses>`_. The numerical factor

    .. math::
        \\frac{\\mathrm{g}}{\\mathrm{cm}^3} =
        \\frac{10^{-15}}{1.67377585}
        \\frac{m_N}{\\mathrm{fm}^3}

    Args:
        rho (float): density :math:`\\rho` [g cm :sup:`-3`]

    Returns:
        float: rho :math:`\\rho` [fm :sup:`-3`]

    See also:
        :func:`fm3togcm3`
    """
    rho = np.asarray(rho)
    return rho/1.67377585e15

RHOSAT = gcm3tofm3(3.e14)
RHOND  = gcm3tofm3(4.e11)

if __name__ == '__main__':
    print("# Unit conversion: 1/fm^3 into g/cm^3")
    print("# 1/fm^3 = 1.67 * 10^15 g/cm^3")
    print("")
    for dens in [1., 0.08, 0.05, 0.03]:
        print("{:f} [1/fm^3] ==  {:.2E} [g/cm^3]".format(dens,fm3togcm3(dens)))

    # print("{0:.2E}".format(fm3togcm3(.08)))
    # print("{0:.2E}".format(fm3togcm3(.05)))
    # print("{0:.2E}".format(fm3togcm3(.03)))

    print("\n\n")
    print("# Unit conversion: g/cm^3 into 1/fm^3.")
    print("# g/cm^3 = 0.59 * 10^-15 1/fm^3")
    print("# Nuclear saturation: {0:.2E}".format(gcm3tofm3(3.e14)))
    print("# Neutron drip: {0:.2E}".format(gcm3tofm3(4.e11)))
    print("{:.6f}".format(gcm3tofm3(4.e11)))
