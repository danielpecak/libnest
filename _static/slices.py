#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Info: who, where, when
# Author: Daniel Pęcak <daniel.pecak@pw.edu.pl>
# Warsaw Technical University, Université Libre de Bruxelles
# On leave: Institute of Physics, Polish Academy of Sciences, Warsaw
# August 2023, Geneva
# =========== Description
# Plotting the cross-sections for different files to check by eye if
# everything is correct and according to the plan.
# -------------------
import numpy as np
import matplotlib.pyplot as plt
import sys
from wdata.io import WData, Var

filenames = ["no_Coulomb_long_06.wtxt"]
filenames = ["no_Coulomb_long_02.wtxt", "no_Coulomb_long_06.wtxt",
"no_Coulomb_long_10.wtxt", "no_Coulomb_long_15.wtxt", "no_Coulomb_long_19.wtxt",
"no_Coulomb_long_24.wtxt", "no_Coulomb_long_27.wtxt",
"no_Coulomb_long_30.wtxt", "no_Coulomb_long_34.wtxt"]


# PATH2   ='/media/data/supercomputing/meff/crust-interaction/wbox/'
PATH2   ="/home/pecak/sshfs/dwarf-scratch/marek/lumi/no_coulomb/init/"

print("# prefix\tN\tP")
for filename in filenames:
    prefix="n=0.0"+filename[16:-5]
    # WBSK data
    data    = WData.load(PATH2+filename, check_data=False)
    [nx,ny,nz] = [data.xyz[i].size for i in range(3)]
    datax = data.xyz[0].reshape(-1,)
    datay = data.xyz[1].reshape(-1,)
    dataz = data.xyz[2].reshape(-1,)

    data.variables[0].filename=PATH2+filename[:-5]+'_'+data.variables[0].name+'.wdat'
    data.variables[1].filename=PATH2+filename[:-5]+'_'+data.variables[1].name+'.wdat'

    dx = datax[1]-datax[0]
    dy = datay[1]-datay[0]
    dz = dataz[1]-dataz[0]
    dv = dx*dy*dz
    # times = np.arange(len(data.t))
    t0=-1 # last time step
    # t0=0 # first time step

    rhon = data.rho_n[t0]
    rhop = data.rho_p[t0]
    # An = data.A_n[t0]
    # Ap = data.A_p[t0]
    Bn = data.B_n[t0]
    Bp = data.B_p[t0]
    Vn = data.V_n[t0]
    Vp = data.V_p[t0]
    Vc = data.V_coulomb[t0]
    Vextn = data.V_ext_n[t0]
    Vextp = data.V_ext_p[t0]
    nun = data.nu_n[t0]
    nup = data.nu_p[t0]
    taun = data.tau_n[t0]
    taup = data.tau_p[t0]
    deltan = data.delta_n[t0]
    deltap = data.delta_p[t0]
    numn = rhon.sum()*dv
    nump = rhop.sum()*dv
    print("{}\t{:10.2f}\t{:.1f}".format(prefix, numn, nump))

    # First create some toy data:
    # x = np.linspace(0, 2*np.pi, 400)
    # y = np.sin(x**2)

    # Create just a figure and only one subplot
    # fig, ax = plt.subplots()
    # ax.plot(x, y)
    # ax.set_title('Simple plot')

    # Create two subplots and unpack the output array immediately
    f, axes = plt.subplots(7, 2, sharex=True)
    ax = axes.flatten()
    f.suptitle('Density: {}'.format(prefix))
    ax[0].set_title('neutronss')
    ax[1].set_title('protons')
    ax[0].set_ylabel(r'$\rho_n$')
    ax[1].set_ylabel(r'$\rho_p$')
    ax[2].set_ylabel(r'$\Delta_n$')
    ax[3].set_ylabel(r'$\Delta_p$')
    ax[4].set_ylabel(r'$\nu_n$')
    ax[5].set_ylabel(r'$\nu_p$')
    ax[6].set_ylabel(r'$B_n$')
    ax[7].set_ylabel(r'$B_p$')
    ax[8].set_ylabel(r'$\tau_n$')
    ax[9].set_ylabel(r'$\tau_p$')
    ax[10].set_ylabel(r'$V_n$')
    ax[11].set_ylabel(r'$V_p$')
    ax[12].set_ylabel(r'$Vext_n$')
    ax[13].set_ylabel(r'$Vext_p$')
    x = dataz
    ax[0].plot(x, rhon[int(nx/2)][int(ny/2)])
    ax[1].plot(x, rhop[int(nx/2)][int(ny/2)])
    ax[2].plot(x, abs(deltan[int(nx/2)][int(ny/2)]))
    ax[3].plot(x, abs(deltap[int(nx/2)][int(ny/2)]))
    ax[4].plot(x, abs(nun[int(nx/2)][int(ny/2)]))
    ax[5].plot(x, abs(nup[int(nx/2)][int(ny/2)]))
    ax[6].plot(x, abs(Bn[int(nx/2)][int(ny/2)]))
    ax[7].plot(x, abs(Bp[int(nx/2)][int(ny/2)]))
    ax[8].plot(x, abs(taun[int(nx/2)][int(ny/2)]))
    ax[9].plot(x, abs(taup[int(nx/2)][int(ny/2)]))
    ax[10].plot(x, abs(Vn[int(nx/2)][int(ny/2)]))
    ax[11].plot(x, abs(Vp[int(nx/2)][int(ny/2)]))
    ax[12].plot(x, abs(Vextn[int(nx/2)][int(ny/2)]))
    ax[13].plot(x, abs(Vextp[int(nx/2)][int(ny/2)]))


    # plt.show()
    plt.savefig('start_{}.png'.format(prefix))
