# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:57:45 2022

@author: aleks
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import libnest.definitions
from libnest.real_data_plots import TXT_PATH, TXT_PATH_UNIFORM, file_check, cross_section_distance, files_set_particles, files_set_type, pairing_field

# TXT_PATH = "C:\\Users\\aleks\\OneDrive\\Dokumenty\\libnest\\txt\\"
# TXT_PATH_UNIFORM = "C:\\Users\\aleks\\OneDrive\\Dokumenty\\libnest\\uniform-txt\\"
#depends on the user (path did not work in main)

# ================================
#      Temperature variation
# ================================
def plot_temperature_delta(particles_nr):
    """
    Returns four graphs illustrating the cross section of the pairing field,
    :math:`\\Delta` [MeV], and temperature [MeV/kB] (where kB is the Boltzmann
    constant). The data is taken from a folder in the directory, the path provided
    earlier (and meant to be defined by the user). The folder contains data files
    for varying amount of particles, so the function accepts the number of particles
    to be plotted as an argument.

    The first graph contains :math:`\\Delta` of all temperatures, the second is
    a closer look at the lower temperatures (first files in order), the third
    represents higher temperatures, and the last graph provides a look at the highest
    temperatures (for which ;math:`\\Delta` forms no obvious pattern).

    Args:
        particles_nr (string): choice of files with a specified number of particles

    Returns:
        None

    See also:
        :func:`cross_section_distance()`
        :func:`pairing_field()`
        :func:`files_set_type(data_type, filenames)`
        :func:`files_set_particles(filenames)`
    """

#MAIN POINTS
# 216 nr of particles:
# - magnitude of delta drops with temperature
# - temperatures 0 and 2 MeV/kB seem to give the exact same pattern
# - there is a huge gap in order of magnitude of delta between 20 and 24 MeV/kB
# - for larger temperatures delta drops quickly (but temperatures of 0-20 MeV/kB
#     were spaced by 2 MeV/kB, while 20-60 by 4 (causing a larger optical difference))
# - delta for temperatures 24 MeV/kB  and larger forms distinct central peaks
#     (24, 28 also form side troughs)
# - delta for temperatures of 52 MeV/kB and above forms no clear pattern

# larger nr of particles:
# - difference between large and low temperatures visible, but not as distinct
#     as with 216 particles: no pairing field close to zero, all have similar shapes
# - magnitude of delta drops with temperature
# - up to around 20 - 30 MeV/kB delta is very similar for all temperatures

    filenames = files_set_type('delta', files_set_particles(particles_nr, TXT_PATH))
    # print(filenames)
    path_filenames = [TXT_PATH + x for x in filenames]
    # print(new_filenames)

    #plotting all temperatures
    plt.figure(figsize=(10,8))
    plt.title(r"Pairing field $\Delta$ varying with temperature", fontsize=15)
    plt.xlabel(r"$ r\: [fm]$", fontsize=10)
    plt.ylabel(r"$\Delta$ [MeV]", fontsize=10)
    plt.xticks(fontsize=10)

    for file in path_filenames:

        if file_check(file):

            data = np.genfromtxt(file, delimiter=' ', comments='#')
            data = data[~np.isnan(data).any(axis=1)]
            data = data[data[:, -1] != 0]
            #data[:,0] - x
            #data[:,1] - y
            #data[:,2] - delta_real
            #data[:,3] - delta_imaginary

            r = cross_section_distance(data[:,0], data[:,1], 180)
            delta, arg = pairing_field(data[:,2], data[:,3])

            plt.scatter(r, delta, 0.5, label=file[-16:-12]+' MeV/kB')

        else:
            sys.exit('# ERROR: Cannot access file')
    plt.legend(title='# particles: '+particles_nr, loc='upper right', ncol=1, markerscale=7)
    plt.show()

    #plotting the top
    plt.figure(figsize=(10,8))
    plt.title(r"Pairing field $\Delta$ varying with temperature - low temperatures", fontsize=15)
    plt.xlabel(r"$ r\: [fm]$", fontsize=10)
    plt.ylabel(r"$\Delta$ [MeV]", fontsize=10)
    # plt.ylim([0.33, 0.35]) #useful for 216 particles
    plt.xticks(fontsize=10)

    top_filenames = path_filenames[0:13]
    for file in top_filenames:

        if file_check(file):

            data = np.genfromtxt(file, delimiter=' ', comments='#')
            data = data[~np.isnan(data).any(axis=1)]
            data = data[data[:, -1] != 0]
            #data[:,0] - x
            #data[:,1] - y
            #data[:,2] - delta_real
            #data[:,3] - delta_imaginary

            r = cross_section_distance(data[:,0], data[:,1], 180)
            delta, arg = pairing_field(data[:,2], data[:,3])

            plt.scatter(r, delta, 0.5, label=file[-16:-12]+' MeV/kB')

        else:
            sys.exit('# ERROR: Cannot access file')
    plt.legend(title='# particles: '+particles_nr, loc='upper right', ncol=1, markerscale=7)
    plt.show()

    #plotting the bottom - can use log scale for lower numbers of particles
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)
    # ax.set_yscale('log') #useful for 216 particles
    plt.title(r"Pairing field $\Delta$ varying with temperature - high temperatures", fontsize=15)
    plt.xlabel(r"$ r\: [fm]$", fontsize=10)
    plt.ylabel(r"$\Delta$ [MeV]", fontsize=10)
    # plt.ylim([1e-15, 5e-3]) # useful for 216 particles
    plt.xticks(fontsize=10)

    bottom_filenames = path_filenames[11:21]
    for file in bottom_filenames:

        if file_check(file):

            data = np.genfromtxt(file, delimiter=' ', comments='#')
            data = data[~np.isnan(data).any(axis=1)]
            data = data[data[:, -1] != 0]
            #data[:,0] - x
            #data[:,1] - y
            #data[:,2] - delta_real
            #data[:,3] - delta_imaginary


            r = cross_section_distance(data[:,0], data[:,1], 180)
            delta, arg = pairing_field(data[:,2], data[:,3])

            plt.scatter(r, delta, 0.5, label=file[-16:-12]+' MeV/kB')

        else:
            sys.exit('# ERROR: Cannot access file')

    plt.legend(title='# particles: '+particles_nr, loc='upper right', ncol=1, markerscale=7)
    plt.show()

    #plotting the highest temperatures (see near zero delta for 216 particles)
    plt.figure(figsize=(10,8))
    plt.title(r"Pairing field $\Delta$ varying with temperature -"
              " highest temperatures", fontsize=15)
    plt.xlabel(r"$ r\: [fm]$", fontsize=10)
    plt.ylabel(r"$\Delta$ [MeV]", fontsize=10)
    plt.xticks(fontsize=10)

    zero_filenames = path_filenames[18:21]
    for file in zero_filenames:

        if file_check(file):

            data = np.genfromtxt(file, delimiter=' ', comments='#')
            data = data[~np.isnan(data).any(axis=1)]
            data = data[data[:, -1] != 0]
            #data[:,0] - x
            #data[:,1] - y
            #data[:,2] - delta_real
            #data[:,3] - delta_imaginary

            r = cross_section_distance(data[:,0], data[:,1], 180)
            delta, arg = pairing_field(data[:,2], data[:,3])

            plt.scatter(r, delta, 0.5, label=file[-16:-12]+' MeV/kB')

        else:
            sys.exit('# ERROR: Cannot access file')
    plt.legend(title='# particles: '+particles_nr, loc='upper right', ncol=1, markerscale=7)
    plt.show()

def plot_max_delta_temperature(particles_nr):
    """
    Plots the maximum value of pairing field :math:`\\Delta` [MeV] against temperature
    T [MeV/k:sub:`B`]. The function creates an array of temperatures (taken from filenames)
    and parses through the files to find the maximum pairing field :math:`\\Delta_{max}`
    in the range 40-60 fm (which is assumed to be the flattest part of the curve).

    Args:
        particles_nr (string): choice of files with a specified number of particles

    Returns
        None
    """
    filenames = files_set_type('delta', files_set_particles(particles_nr, TXT_PATH))
    path_filenames = [TXT_PATH + x for x in filenames]
    delta_max = []
    temperature = []

    for file in path_filenames:
        if file_check(file):
            data = np.genfromtxt(file, delimiter=' ', comments='#')
            data = data[~np.isnan(data).any(axis=1)]
            data = data[data[:, -1] != 0]
            #data[:,0] - x
            #data[:,1] - y
            #data[:,2] - delta_real
            #data[:,3] - delta_imaginary

            r = cross_section_distance(data[:,0], data[:,1], 180)
            i = np.where(np.logical_and(r>=40, r<=60))
            delta, arg = pairing_field(data[i,2], data[i,3])
            delta_max.append(np.max(delta))
            temperature.append(float(file[-16:-12]))
        else:
            sys.exit('# ERROR: Cannot access file')

    plt.figure()
    plt.title(r"Pairing field $\Delta$ for "+particles_nr+" particles", fontsize=15)
    plt.xlabel(r"$ T \: [MeV/k_B]$", fontsize=10)
    plt.ylabel(r"$\Delta_{max}$ [MeV]", fontsize=10)
    plt.plot(temperature, delta_max)
    # plt.legend(title='# particles: '+particles_nr, loc='upper right', ncol=1, markerscale=7)
    plt.show()


def plot_max_delta_temperature_uniform(particles_nr):
    """
    Plots the maximum value of pairing field :math:`\\Delta` [MeV] against temperature
    T [MeV/k:sub:`B`]. The function creates an array of temperatures (taken from filenames)
    and parses through the files to find the maximum pairing field :math:`\\Delta_{max}`
    in the range 40-60 fm (which is assumed to be the flattest part of the curve).

    Args:
        particles_nr (string): choice of files with a specified number of particles

    Returns
        None
    """
    filenames = files_set_type('delta', files_set_particles(particles_nr, TXT_PATH))
    path_filenames = [TXT_PATH + x for x in filenames]
    delta_max = []
    temperature = []

    for file in path_filenames:
        if file_check(file):
            data = np.genfromtxt(file, delimiter=' ', comments='#')
            data = data[~np.isnan(data).any(axis=1)]
            data = data[data[:, -1] != 0]
            #data[:,0] - x
            #data[:,1] - y
            #data[:,2] - delta_real
            #data[:,3] - delta_imaginary

            r = cross_section_distance(data[:,0], data[:,1], 180)
            # i = np.where(np.logical_and(r>=40, r<=60))
            i = np.where(np.logical_and(r>=0, r<=10))
            delta, arg = pairing_field(data[i,2], data[i,3])
            delta_max.append(np.max(delta))
            temperature.append(float(file[-16:-12]))
        else:
            sys.exit('# ERROR: Cannot access file')
    
    filenames_uniform = files_set_type('delta', files_set_particles(particles_nr, TXT_PATH_UNIFORM))
    path_filenames_uniform = [TXT_PATH_UNIFORM + x for x in filenames_uniform]
    delta_max_uniform = []
    temperature_uniform = []

    for file in path_filenames_uniform:
        if file_check(file):
            data = np.genfromtxt(file, delimiter=' ', comments='#')
            data = data[~np.isnan(data).any(axis=1)]
            data = data[data[:, -1] != 0]
            #data[:,0] - x
            #data[:,1] - y
            #data[:,2] - delta_real
            #data[:,3] - delta_imaginary

            r = cross_section_distance(data[:,0], data[:,1], 180)
            # i = np.where(np.logical_and(r>=40, r<=60))
            i = np.where(np.logical_and(r>=0, r<=10))
            delta, arg = pairing_field(data[i,2], data[i,3])
            delta_max_uniform.append(np.max(delta))
            temperature_uniform.append(float(file[-16:-12]))
        else:
            sys.exit('# ERROR: Cannot access file')

    # temperature_critical = delta_max/1.764
    delta_max = np.array(delta_max)
    # temperature_critical = np.mean(delta_max[0:4]/1.764)
    temperature_critical = delta_max_uniform[0]/1.764
    print(temperature_critical)

    delta_theoretical_t_critical = 3.06*temperature_critical*np.sqrt(1-temperature/temperature_critical)
    #np.where((temperature>temperature_critical), 3.06*temperature_critical*np.sqrt(1-temperature/temperature_critical), delta_theoretical_t_critical)
    x = np.array(temperature)/float(temperature_critical)
    x_uniform = np.array(temperature_uniform)/float(temperature_critical)

    plt.figure()
    plt.xlim(0., 1.1)
    plt.title(r"Pairing field $\Delta$ for "+particles_nr+" particles", fontsize=15)
    plt.xlabel(r"$ T/T_{crit}$", fontsize=10)
    plt.ylabel(r"$\Delta_{max}$ [MeV]", fontsize=10)
    plt.plot(x, delta_max, '-o', label="vortex")
    plt.plot(x_uniform, delta_max_uniform, '-o', label="uniform")
    plt.plot(x[4:], delta_theoretical_t_critical[4:], label=r"delta for $T \to T_{crit}$", linestyle='dashed')
    plt.legend(ncol=1)
    plt.show()
    

    
if __name__ == '__main__':
    pass
