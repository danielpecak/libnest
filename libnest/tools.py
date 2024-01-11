#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Info: who, where, when
# @author Daniel Pęcak <Daniel.Pecak@pw.edu.pl>
# Warsaw Technical University
# On leave: Institute of Physics, Polish Academy of Sciences, Warsaw
# October 2022, Brussels
"""
tools.py
========
Module for analyzing the wdata for the inner crust.
That includes slices through 3D and 2D data and genating
2D and 1D data of reduced dimensionality.

"""
import numpy as np
import math
from libnest import units
from libnest.units import DENSEPSILON
from libnest.definitions import rho2kf
from libnest.bsk import eF_n

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
        :func:`kf2rho`
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

    Args:
        density numpy array

    Returns:
        2 dimensional numpy array of center of mass coordinates
        at subsequent time steps [[x0,y0,z0], [x1,y1,z1],...]
    """
    [nx, ny, nz] = [density.shape[1], density.shape[2], density.shape[3]]
    position=[]
    for t in range(density.shape[0]):
        mass=density[t]
        total_mass=np.sum(mass)
        x=0
        y=0
        z=0
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x+=i*mass[i, j, k]
                    y+=j*mass[i, j, k]
                    z+=k*mass[i, j, k]
        position.append([x/total_mass, y/total_mass, z/total_mass])
    return np.array(position)

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
    [nx, ny, nz] = [density.shape[1], density.shape[2], density.shape[3]]
    energy=[]
    for t in range(density.shape[0]):
        e=0
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    e+=np.absolute(delta[t][i, j, k])**2/(eF_n(rho2kf(density[t][i, j ,k]))*density[t][i, j, k])
        energy.append(3./8.*e)
    return np.array(energy)


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
    [nx, ny, nz] = [density_n.shape[1], density_n.shape[2], density_n.shape[3]]
    energy=[]
    density=density_n+density_p
    for t in range(density.shape[0]):
        e=0
        for i in range(nx):
            for m in range(ny):
                for k in range(nz):
                    e+=(j[t][0][i, 0, 0]**2+j[t][1][0, m, 0]**2+j[t][2][0, 0, k]**2)/(density[t][i, m, k])

        energy.append(e*HBARC/(2*(MN*particleN(density_n)[t]+MP*particleN(density_p)[t]))/(particleN(density)[t]))

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
    particles=[]
    for t in range(density.shape[0]):
        particles.append(np.sum(density[t]))

    return np.array(particles)


if __name__ == '__main__':
    pass
