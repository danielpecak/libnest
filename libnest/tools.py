#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Info: who, where, when
# @author Daniel Pęcak <Daniel.Pecak@pw.edu.pl>
# Warsaw Technical University
# On leave: Institute of Physics, Polish Academy of Sciences, Warsaw
# October 2022, Brussels
"""
Module: tools
=============
List of functions
-----------------
"""
import numpy as np
import math
from libnest import units
from libnest.units import DENSEPSILON
from libnest.definitions import rho2kf, eF_n

def threeSlice(variable):
    """
    Returns 1, 2 or 3 dimensional numpy array with slices through given variable.

    It works for 1D, 2D, 3D.

    By default it cuts through the center.


    Args:
        numpy array

    Returns:
        n dimensional numpy array with slices

    See also:
        :func:`.kf2rho`
    """
    if(len(variable.shape)==3):
        [nx, ny, nz] = [variable.shape[i] for i in range(3)]
        return np.asarray([variable[:, int(ny/2),int(nz/2)],
                         variable[int(nx/2), :, int(nz/2)],
                         variable[int(nx/2),int(ny/2), :]], dtype="object")
    elif(len(variable.shape)==2):
        [nx,ny] = [variable.shape[i] for i in range(2)]
        return np.asarray(variable[:,int(ny/2)], variable[int(nx/2),:], dtype="object")
    elif(len(variable.shape)==1):
        nx=variable.shape[0]
        return np.asarray([variable], dtype="object")
    return 0




def centerOfMass(density):
    """
    Returns numpy array of center of masses coordinates.

    It works for 1D, 2D, 3D.

    .. note::
        The center of mass unit here is unitless (units of grid)

    Args:
        density(float): :math:`\\rho_n` [fm :sup:`-3`]

    Returns:
        2 dimensional numpy array of center of mass coordinates
        at subsequent time steps [[x0,y0,z0], [x1,y1,z1],...]
    """
    i, j, k = np.indices(density.shape[1:])
    total_mass = np.sum(density, axis=(1, 2, 3))
    x = np.sum(i * density, axis=(1, 2, 3))
    y = np.sum(j * density, axis=(1, 2, 3))
    z = np.sum(k * density, axis=(1, 2, 3))
    position = np.column_stack((x / total_mass, y / total_mass, z / total_mass))
    return position

def condensationEnergy(density, delta,dV=1.0):
    """
    Calculates the condensation energy by calculating proper function of system density :math:`\\rho`, and pairing field :math:`\\Delta`. The Fermi energy :math:`\\epsilon^*_{F}` is calculated under the assumption that the effective mass is taken into account.
    It works for 1D, 2D, 3D.
    Returns numpy array of values of condensation energy.
    1 dimensional numpy array of condensation energy


    .. math::

	   E_{\\mathrm{cond}} = \\int_V \\frac{3}{8} \\frac{|\\Delta(\\bm r)|^2}{\\epsilon^*_{F}(\\bm r)}\\rho(\\bm r)d{\\bm r}


    Args:
        density (float): :math:`\\rho` [fm :sup:`-3`]
        delta (float): :math:`\\Delta` [MeV]
        dV (float): infinitesimal element of volume (by default set to 1.0)

    Returns:
        float: condensation energy :math:`E_{\\mathrm{cond}}` at subsequent time steps [e0, e1,...]

    """
    e=np.sum(3./8.*np.absolute(delta)**2/(eF_n(rho2kf(density))*(density)+10**-12), axis=(1, 2, 3))
    return e


def flowEnergy(j, density_n, density_p):
    """
    Returns numpy array of values of flow energy.


    It works for 1D, 2D, 3D.

    Formula: \integral hbar*c*j^2/(2mc^2 density) dr

    Args:
        current numpy array, density of neutrons numpy array, density of protons numpy array

    Returns:
        1 dimensional numpy array of flow energy
        at subsequent time steps [e0, e1,...]

    """
    density=density_n+density_p
    i, m, k = np.indices(density.shape[1:])
    e=np.sum((j[:,0]**2+j[:, 1]**2+j[:, 2]**2)/(density), axis=(1,2,3))
    energy=e*HBARC/(2*(MN*particleN(density_n)+MP*particleN(density_p))/(particleN(density)))
    return np.array(energy)



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
    return np.sum(density, axis=(1, 2, 3))


if __name__ == '__main__':
    pass
