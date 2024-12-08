#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Description
#
#================================
"""
nucleus.py
======
The module that contains toolkit for dealing for nuclei, including
nuclei immersed in the neutron matter.

==========

"""
import sys
import numpy as np
# from libnest import units
# from libnest.units import HBARC, DENSEPSILON, NUMZERO
# from libnest.units import MN, MP, HBAR2M_n, HBAR2M_p
# from libnest.definitions import rho2kf, rhoEta, rho2tau

# ================================
#       ???
# ================================
def cm(data,t=-1):
    """
    Center of mass
    """
    if(type(t)==int):
        t = np.array([t])
    # t = np.array(t)
    res = []
    [nx,ny,nz] = [data.xyz[i].size for i in range(3)]
    datax = data.xyz[0].reshape(-1,)
    datay = data.xyz[1].reshape(-1,)
    dataz = data.xyz[2].reshape(-1,)
    dx = datax[1] - datax[0]
    dy = datay[1] - datay[0]
    dz = dataz[1] - dataz[0]
    dv = dx*dy*dz

    V  = dv*nx*ny*nz
    rho_back = data.rho_n[0][int(nx/2)][int(ny/2)][nz-1]
    # print(rho_back)
    for t0 in t:
        q = [0.0, 0.0, 0.0]
        qx = 0.0
        qy = 0.0
        qz = 0.0
        # TODO DEALING WITH background density
        rhon = data.rho_n[t0] - rho_back
        rhop = data.rho_p[t0]
        # rho = rhop + rhon
        rho = rhop
        nn = 0*rhon.sum() # NOTE missing dv, it will cancel out
        pp = rhop.sum() # NOTE missing dv, it will cancel out
        na = nn + pp    # NOTE missing dv, it will cancel out
        # print(nn*dv,np*dv)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    # TODO include X, Y, Z directions
                    qx = qx + rho[i][j][k]*datax[i]
                    qy = qy + rho[i][j][k]*datay[j]
                    qz = qz + rho[i][j][k]*dataz[k]
        q[0] = qx/na
        q[1] = qy/na
        q[2] = qz/na
        res.append(q)
    if len(res)==1:
        return(res[0])
    else:
        return(res)
    # return(0.0)

# ================================
#       Deformation
# ================================
def q20(data,t=-1,rcm=np.array([[0,0,0]])):
    """
    Quadrupole moment

    .. math::

	   Q_{20} = \\int

    Args:
        rho_n (float): neutron density :math:`\\rho_n` [fm :sup:`-3`]; sum of both spin components

    Returns:
        float: pairing field for neutron matter :math:`\\Delta_{\\mathrm{NeuM}}` [MeV]

    See also:
        :func:`symmetric_pairing_field`
        :func:`neutron_ref_pairing_field`
        :func:`proton_ref_pairing_field`
    """
    # TODO DODUMENTATION, FORMULA
    # data # wdata
    # t0=-1 # if not specified the last timestep
    if(type(t)==int):
        t = [t]
    res = []
    # print(rcm)
    rcm = np.array(rcm)
    # print(rcm.shape)
    [nx,ny,nz] = [data.xyz[i].size for i in range(3)]
    dataxt = data.xyz[0].reshape(-1,)
    datayt = data.xyz[1].reshape(-1,)
    datazt = data.xyz[2].reshape(-1,)
    dx = dataxt[1] - dataxt[0]
    dy = datayt[1] - datayt[0]
    dz = datazt[1] - datazt[0]
    dv = dx*dy*dz
    V  = dv*nx*ny*nz
    for t0 in range(len(t)):
        print("# Timestep: {}".format(t[t0]))
        # print(rcm[t0])
        q = [0.0, 0.0, 0.0]
        qx = 0.0
        qy = 0.0
        qz = 0.0
        # TODO DEALING WITH background density
        rhon = data.rho_p[t[t0]]
        datax = dataxt - rcm[t0][0]
        datay = datayt - rcm[t0][1]
        dataz = datazt - rcm[t0][2]
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    qx = qx + rhon[i][j][k]*datax[i]*datax[i]
                    qy = qy + rhon[i][j][k]*datay[j]*datay[j]
                    qz = qz + rhon[i][j][k]*dataz[k]*dataz[k]
        q[0] = (2*qx - qy - qz)/(qx+qy+qz)*np.sqrt(np.pi/5)
        q[1] = (2*qy - qx - qz)/(qx+qy+qz)*np.sqrt(np.pi/5)
        q[2] = (2*qz - qy - qx)/(qx+qy+qz)*np.sqrt(np.pi/5)
        # q[0] = (2*qx - qy - qz)*dv#/(qx+qy+qz)
        # q[1] = (2*qy - qx - qz)*dv#/(qx+qy+qz)
        # q[2] = (2*qz - qy - qx)*dv#/(qx+qy+qz)
        res.append(q)
    if len(res)==1:
        return(res[0])
    else:
        return(res)



if __name__ == '__main__':
    pass
