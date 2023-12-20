#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Info: who, where, when
# @author Daniel Pęcak <Daniel.Pecak@pw.edu.pl>
# Warsaw Technical University
# On leave: Institute of Physics, Polish Academy of Sciences, Warsaw
# October 2022, Brussels
"""
tools.py
==============
"""
import numpy as np
import math
from libnest import units
from libnest.units import DENSEPSILON


def threeSlice(file):
    """
    Returns 1, 2 or 3 numpy arrays with slices.

    By default it cuts through the center.


    Args:
        file (WData):

    Returns:
        float: TODO

    See also:
        :func:`kf2rho`
    """
    #TODO
    return 0



def centerOfMass(density):
    """
    Returns numpy array of center of masses coordinates.

    It works for 1D, 2D, 3D.

    Args:
        density numpy array

    Returns:
        2 dimensional numpy array of center of mass coordinates
        at subsequent time steps [[x0,y0,z0], [x1,y1,z1],...]

    """
    # TODO: MM
    return 0

def condensationEnergy(density, delta):
    """
    Returns numpy array of values of condensation energy.

    It works for 1D, 2D, 3D.

    Formula: integral 3/8 |Delta(r)|^2/EFermi(r)density(r)dr

    Args:
        density numpy array, delta numpy array

    Returns:
        1 dimensional numpy array of condensation energy
        at subsequent time steps [e0, e1,...]

    """
    # TODO: MM
    return 0


def flowEnergy(j, density):
    """
    Returns numpy array of values of flow energy.

    It works for 1D, 2D, 3D.

    Formula: \integral hbar*c*j^2/(2mc^2 density) dr

    Args:
        density numpy array, delta numpy array

    Returns:
        1 dimensional numpy array of condensation energy
        at subsequent time steps [e0, e1,...]

    """
    # TODO: MM
    return 0



def particleN(density):
    """
    Returns number of particles.

    It works for 1D, 2D, 3D.
    For static nucleus it is constant, however there are some system
    scenarios where the number of particles might change.

    Args:
        density numpy array

    Returns:
        1 dimensional numpy array of number of particles
        at subsequent time steps [N0, N1,...]

    """
    # TODO: MM
    return 0


if __name__ == '__main__':
    pass
