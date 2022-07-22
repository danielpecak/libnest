#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Info: who, where, when
# @author Daniel Pęcak <Daniel.Pecak@pw.edu.pl>
# Warsaw Technical University, Université Libre de Bruxelles
# On leave: Institute of Physics, Polish Academy of Sciences, Warsaw
# March 2022, Brussels
"""
units.py
========
Module for converting the density units for nuclear matter in the ranges typical
for neutron stars.

Numerical
---------
Constants for dealing with extreme numerics:
 * ``NUMZERO`` = 1e-12 (the lowest numerical value allowed)
 * ``NUMINF``  = 1e30  (the largest numerical value allowed)
 * ``DENSEPSILON`` = 1e-12

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
MN        939.5654205 [MeV]               neutron mass
MP        938.2720882 [MeV]               proton  mass
HBAR2M_n  20.72124837 [MeV*fm :math:`^2`] 0.5*hbar^2/MN
HBAR2M_p  20.74981092 [MeV*fm :math:`^2`] 0.5*hbar^2/MP
========  =========== =================== ===========

"""
import numpy as np

DENSEPSILON = 1e-12
NUMZERO = 1e-12
NUMINF  = 1e30


HBARC=197.3269804  # \hbar c [MeV fm]
VUNIT=197./940*100 # alpha constant over the neutron mass (100 to have percentage)
hbar22M0   =20.72  # neutron bare mass
MN   =939.56542052 # neutron mass [MeV]
MP   =938.27208816 # proton  mass [MeV]
HBAR2M_n = 20.721248369006936 # [MeV*fm<sup>2</sup>] 0.5*hbar^2/neutron mass
HBAR2M_p = 20.749810921501915 # [MeV*fm<sup>2</sup>] 0.5*hbar^2/proton mass


def KtoMev(temp):
    """Function converts temperature units from Kelvins to MeVs.


    .. todo::
        formula

    Args:
        temp (float): temperature :math:`T` [K]

    Returns:
        float: temperature :math:`T` [MeV]

    See also:
        :func:`MeVtoK`
    """
    return temp/11604525006.1598

def MeVtoK(temp):
    """Function converts temperature units from MeVs to Kelvins.


    .. todo::
        formula

    Args:
        temp (float): temperature :math:`T` [MeV]

    Returns:
        float: temperature :math:`T` [K]

    See also:
        :func:`KtoMev`
    """

    return temp*11604525006.1598

def fm3togcm3(rho):
    """Function converts units 1/fm^3 into g/cm^3.

    .. todo::
        formula

    Args:
        temp (rho): density :math:`\\rho` [1/fm^3]

    Returns:
        float: rho :math:`\\rho` [g/cm^3]

    See also:
        :func:`gcm3tofm3`
    """
    return rho*1.67e15

def gcm3tofm3(rho):
    """Function converts units g/cm^3 into 1/fm^3.

    .. todo::
        formula

    Args:
        temp (rho): density :math:`\\rho` [g/cm^3]

    Returns:
        float: rho :math:`\\rho` [1/fm^3]

    See also:
        :func:`fm3togcm3`
    """
    return rho/1.67e15

if __name__ == '__main__':
    print("# Unit conversion: 1/fm^3 into g/cm^3")
    print("# 1/fm^3 = 1.67 * 10^15 g/cm^3")
    print("{0:.2E}".format(fm3togcm3(1.)))
    print("{0:.2E}".format(fm3togcm3(.08)))
    print("{0:.2E}".format(fm3togcm3(.05)))
    print("{0:.2E}".format(fm3togcm3(.03)))

    print("\n\n")
    print("# Unit conversion: g/cm^3 into 1/fm^3.")
    print("# g/cm^3 = 0.59 * 10^-15 1/fm^3")
    print("# Nuclear saturation: {0:.2E}".format(gcm3tofm3(3.e14)))
    print("# Neutron drip: {0:.2E}".format(gcm3tofm3(4.e11)))
    print("{:.6f}".format(gcm3tofm3(4.e11)))

rhoSAT = gcm3tofm3(3.e14)
rhoND  = gcm3tofm3(4.e11)
