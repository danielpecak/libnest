#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Info: who, where, when
# Author: Daniel Pęcak <daniel.pecak@pw.edu.pl>
# Warsaw Technical University
# On leave: Institute of Physics, Polish Academy of Sciences, Warsaw
# January 2024, Warsaw
# =========== Description
# Here we test TOOLS module: centre of mass of the neutron density in time,
# saved to crossectionNEW.png.
# =========== Usage example (requires py-WDATA: pip install wdata)
# $ ./tools_example.py path/to/run.wtxt
import sys
import numpy as np
import matplotlib.pyplot as plt
from wdata.io import WData

from libnest import tools
# print( "test")
# print(dir(tools))

if len(sys.argv) != 2:
    sys.exit("Usage: tools_example.py FILE.wtxt")
file = sys.argv[1]
data = WData.load(file, check_data=False)
[nx,ny,nz] = [data.xyz[i].size for i in range(3)]
datax = data.xyz[0].reshape(-1,)  # Array of x coordinates
datat = data.t
# print(datat)
rhon=data.rho_n#[490:]
datat   = datat#[490:]
cms = tools.centerOfMass(rhon)
cms=np.transpose(cms)


t = -1 # Take the last 'timestep'. In our case the larger number of 'time',
       # the more accurate plot

########
## Plotting densities
########
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 16))
ax = ax.flatten()

labels = [r'$\rho_q$', r'$\rho_q$', r'$u_q$', r'$u_q$', r'$|\Delta_q|$', r'$|\Delta_q|$', r'$\nu_q$', r'$\nu_q$', r'$|j_{iq}|$', r'$|j_{iq}|$']
i=0
ax[i].plot(datat, cms[:][2],'-',label=r"$cm_z$")
i+=1
ax[i].plot(datat, cms[:][0],'-',label=r"$cm_x$")
ax[i].plot(datat, cms[:][1],'-',label=r"$cm_y$")
i+=1

# Positioning of labels
for i in range(2):
    ax[i].set_ylabel(labels[i])
    ax[i].legend(loc='upper left')

plt.savefig("crossectionNEW.png")
# plt.show()
plt.close()
