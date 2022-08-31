# -*- coding: utf-8 -*-
"""
Plotting real data
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import libnest.definitions
import libnest.plots

TXT_PATH = "C:\\Users\\aleks\\OneDrive\\Dokumenty\\libnest\\txt\\"
TXT_PATH_UNIFORM = "C:\\Users\\aleks\\OneDrive\\Dokumenty\\libnest\\uniform-txt\\"
#depends on the user (path did not work in main)

# ================================
#          Handling Files
# ================================

def file_check(filename):
    """
    Checks if all files are present in the directory.

    Args:
        filename (string): name of the data set file

    Raises:
        FileNotFoundError

    Returns:
        bool
    """
    try:
        file = open(filename, 'r')
        file.close()
        return True
    except FileNotFoundError:
        print(f"{filename}" + " not found. Check the directory and the file.")
        return False

def files_set_particles(particles_nr, directory_path):
    """
    Returns an array with files containing data sets only for a specified
    number of particles.

    Args:
    particles_nr (string): number of particles of the data set

    Returns
        array: list of filenames
    """
    filenames = []

    for path in os.listdir(directory_path):
        if os.path.isfile(os.path.join(directory_path, path)):
            filenames.append(path)

    p216 = []
    p4160 = []
    p8000 = []
    p13600 = []
    p16640 = []
    p24000 = []

    p216 = [i for i in filenames if '216' in i]
    p4160 = [i for i in filenames if '4160' in i]
    p8000 = [i for i in filenames if '8000' in i]
    p13600 = [i for i in filenames if '13600' in i]
    p16640 = [i for i in filenames if '16640' in i]
    p24000 = [i for i in filenames if '24000' in i]

    if particles_nr == '216':
        return p216
    elif particles_nr == '4160':
        return p4160
    elif particles_nr == '8000':
        return p8000
    elif particles_nr == '13600':
        return p13600
    elif particles_nr == '16640':
        return p16640
    elif particles_nr == '24000':
        return p24000
    else:
        sys.exit('# ERROR: Incorrect value')

def files_set_type(data_type, filenames):
    """
    Takes an array of files and returns an array with files containing data sets of
    a specified kind: rho (matter density), current, delta (pairing field),
    or A (mean field potential)

    Args:
    data_type (string): choice of type of files
    filenames (array): list of filenames to sort through

    Returns
        array: list of filenames of a specified kind
    """
    p_density = []
    p_current = []
    p_delta = []
    p_A = []

    p_density = [i for i in filenames if 'density' in i]
    p_current = [i for i in filenames if 'current' in i]
    p_delta = [i for i in filenames if 'delta' in i]
    p_A = [i for i in filenames if 'A' in i]


    if data_type == 'density':
        return p_density
    elif data_type == 'current':
        return p_current
    elif data_type == 'delta':
        return p_delta
    elif data_type == 'A':
        return p_A
    else:
        sys.exit('# ERROR: Incorrect data type input')

# ================================
#           Coordinates
# ================================

def cross_section_distance(x, y, size):
    """
    Returns the distance to the centre of the 90x90 fm box as cross-section,
    changing it from a coordinate system with origin at one corner of the box.

     Args:
        x (float): x-coordinate of :math:`r` [fm]
        y (float): y=coordinate of :math:`r` [fm`]
        size (float): full size of one axis of the square/box [fm`]

    Returns:
        float: cross section distance :math:`r` [fm]

    """
    r = np.sqrt((x-size/2)**2 + (y-size/2)**2)
    i = np.where(x<size/2)
    r[i] = -r[i]
    return r

def phi(x,y, size): #just in case
    """
    Returns the angle of the line connecting the coordinates to the centre
    of the 90x90 "box", for a coordinate system with origin at a corner of the
    box.

     Args:
        x (float): x-coordinate :math:`r` [fm]
        y (float): y=coordinate :math:`r` [fm`]
        size (float): full size of one axis of the square/box [fm`]

    Returns:
        float: angle :math:`\\phi` [rad]
    """
    n = size/2
    return np.arctan(y-n,x-n)

# ================================
#      Real data definitions
# ================================

def pairing_field(rel, im):
    """
    Calculates the absolute value/modulus and the argument of the field
    potential from its imaginary and real parts.

    Args:
        rel (float): real part of field potential :math:`\\Delta_{rel}` [MeV]
        im (float): imaginary part of field potential :math:`\\Delta_{im}` [MeV]

    Returns:
        float: pairing field :math:`\\Delta` [MeV]
    """
    return np.sqrt(im**2+rel**2), np.arctan(im,rel)

def vector_magnitude(x, y, z):
    """
    Calculates the magnitude of a vector with its vector x, y, and z components.

    Args:
        x (float): :math:`\\vec x` component [fm :sup:`-4`]
        y (float): :math:`\\vec y` component [fm :sup:`-4`]
        z (float): :math:`\\vec z` component [fm :sup:`-4`]

    Returns:
        float: magnitude of a vector [fm :sup:`-4`]

    """
    return np.sqrt(x**2+y**2+z**2)

def jsum2(i_x, i_y, i_z, j_x, j_y, j_z):
    """
    Sum of of currents :math:`\\vec i` and :math:`\\vec j`, squared.

    Args:
        i_x (float): x component of :math:`\\vec i`
        i_y (float): y component of :math:`\\vec i`
        i_z (float): z component of :math:`\\vec i`
        j_x (float): x component of :math:`\\vec j`
        j_y (float): y component of :math:`\\vec j`
        j_z (float): z component of :math:`\\vec j`

    Returns
        float: squared sum of currents, :math:`(\\vec i + \\vec j)^2`
    """
    x = i_x + j_x
    y = i_y + j_y
    z = i_z + j_z
    return x**2 + y**2 + z**2

def jdiff2(i_x, i_y, i_z, j_x, j_y, j_z):
    """
    Difference of currents :math:`\\vec i` and :math:`\\vec j`, squared.

    Args:
        i_x (float): x component of :math:`\\vec i`
        i_y (float): y component of :math:`\\vec i`
        i_z (float): z component of :math:`\\vec i`
        j_x (float): x component of :math:`\\vec j`
        j_y (float): y component of :math:`\\vec j`
        j_z (float): z component of :math:`\\vec j`

    Returns
        float: squared difference of currents :math:`(\\vec j - \\vec j)^2`
    """
    x = i_x - j_x
    y = i_y - j_y
    z = i_z - j_z
    return x**2 + y**2 + z**2

# ================================
#       Plotting functions
# ================================
def plot_density(filename):
    """
    Opens the specified file, checks its validity, and creates a 2D array from
    its data. The first and second columns are x- and y-coordinates. The third
    column is the density :math:`\\rho` [fm :sup:`-3`].

    Cross section distance is calculated and the density is divided by bulk
    density. The ratio is plotted against the cross section distance.

    Args:
        filename (string): name of the data set file

    Returns:
        None

    See also:
        :func:`cross_section_distance`
    """
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - rho_q

        r = cross_section_distance(data[:,0], data[:,1], 90)
        rho_bulk = 0.00590448
        rho_ratio = data[:,2]/rho_bulk

        rho = plt.figure()
        rho.add_subplot(111)
        plt.title("Density vs radius", fontsize=15)
        plt.xlabel(r"$ r \: [fm]$", fontsize=10)
        plt.ylabel(r"$ \rho / \rho_{\infty}$", fontsize=10) #\: {[fm}^{-3}]$
        plt.xticks(fontsize=10)
        plt.scatter(r, rho_ratio, 0.5) #plotting
        #plt.legend()
        plt.show()
    else:
        sys.exit('# ERROR: Cannot access file')

def plot_density_contour(filename):
    """
    Opens the specified file, checks its validity, and creates a 2D array from
    its data. The first and second columns are x- and y-coordinates. The third
    column is the density :math:`\\rho` [fm :sup:`-3`].

    //A grid is created and populated using triangulation ("scipy.interpolate.triangulate").
    A contour plot of density as a vertical cross section is created.

    Args:
        filename (string): name of the data set file

    Returns:
        None
    """
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - rho_q
        # x in one direction, y in another??

        x = data[:,0]
        y = data[:,1]

        xx, yy = np.meshgrid(x, y)
        rho_bulk = 0.00590448
        rho_ratio = data[:,2]/rho_bulk

        #interpolation by triangulation
        zi = scipy.interpolate.griddata((x, y), rho_ratio, (xx, yy), method='linear')
        plt.imshow(zi, vmin=rho_ratio.min(), vmax=rho_ratio.max(), origin='lower',
                   extent=[x.min(), x.max(), y.min(), y.max()])
        plt.colorbar()
        plt.show()

        #interpolation by radial basis function
        # inefficient for larger arrays

        # rbf = scipy.interpolate.Rbf(x, y, rho_ratio, function='linear')
        # zi = rbf(xx, yy)
        # plt.imshow(zi, vmin=rho_ratio.min(), vmax=rho_ratio.max(), origin='lower',
        #    extent=[x.min(), x.max(), y.min(), y.max()])
        # plt.colorbar()
        # plt.show()

        #rho = plt.figure()
        #rho_contour = rho.add_subplot(111)
        # plt.title("Density vs radius", fontsize=15)
        # plt.xlabel(r"$ r \: [fm]$", fontsize=10)
        # plt.ylabel(r"$ r \: [fm]$", fontsize=10) #\: {[fm}^{-3}]$
        # cp = rho_contour.contour(xx, yy, rho_ratio, 1)
        # cp.colotbar()
        # cp.clabel()
        # plt.xticks(fontsize=10)
        # plt.legend()
        # plt.show()

    else:
        sys.exit('# ERROR: Cannot access file')

def plot_pairing_field(filename):
    """
    Opens the specified file, checks its validity, and creates a 2D array from
    its data. The first and second columns are x- and y-coordinates. The third
    column is the real part of the pairing field :math:`\\Delta_{rel}` [MeV] and
    the fourth is the imaginary part of the pairing field :math:`\\Delta_{im}` [MeV].

    Cross section distance and the total pairing field are calculated. The total
    pairing field is divided by the bulk pairing field, and the ratio is plotted
    against the cross section distance.

    Args:
        filename (string): name of the data set file

    Returns:
        None

    See also:
        :func:`cross_section_distance`
        :func:`pairing_field`
    """
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - delta_rel
        #data[:,3] - delta_im

        r = cross_section_distance(data[:,0], data[:,1], 90)
        delta, arg = pairing_field(data[:,2], data[:,3])
        delta_bulk = 1.33394659
        delta_ratio = delta/delta_bulk

        plot = plt.figure()
        plot.add_subplot(111)
        plt.title("Pairing field vs radius", fontsize=15)
        plt.xlabel(r"$ r\: [fm]$", fontsize=10)
        plt.ylabel(r"$\Delta / \Delta_{\infty}$", fontsize=10)
        plt.xticks(fontsize=10)
        plt.scatter(r, delta_ratio, 0.5)
        #plt.legend()
        plt.show()
    else:
        sys.exit('# ERROR: Cannot access file')

def plot_current(filename):
    """
    Opens the specified file, checks its validity, and creates a 2D array from
    its data. The first and second columns are x- and y-coordinates. The third,
    fourth, and fifth columns are :math:`\\vec x`, :math:`\\vec y`, :math:`\\vec z`
    current components [fm :sup:`-4`]

    Cross section distance and the total current vector :math:`\\vec j` are
    calculated. The current is then plotted against the cross section distance.

    Args:
        filename (string): name of the data set file

    Returns:
        None

    See also:
        :func:`cross_section_distance`
        :func:`current`
    """
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        print(data)
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - x current component
        #data[:,3] - y current component
        #data[:,4] - z current component

        r = cross_section_distance(data[:,0], data[:,1], 240)
        print(r)
        j = vector_magnitude(data[:,2], data[:,3], data[:,4])
        print(j)

        #NOTE: treats some files as empty? to double check
        # plots #N960i.3_current.txt but its size is 240, not 90
        # add size to function arguments?

        plot = plt.figure()
        plot.add_subplot(111)
        plt.title("Current vs radius", fontsize=15)
        plt.xlabel(r"$ r\: [fm]$", fontsize=10)
        plt.ylabel("j(r)", fontsize=10)
        plt.xticks(fontsize=10)
        plt.scatter(r, j, 0.5)
        #plt.legend()
        plt.show()
    else:
        sys.exit('# ERROR: Cannot access file')


# ================================
#        All data plots
# ================================
def plot_density_slice(filename):
    """
    Opens the specified file, checks its validity, and creates a 2D array from
    its data. The first and second columns are x- and y-coordinates. The third
    column is the density :math:`\\rho` [fm :sup:`-3`].

    Args:
        filename (string): name of the data set file

    Returns:
        None

    See also:
        :func:`cross_section_distance`
    """
    filename = TXT_PATH + filename
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - rho_q

        r = cross_section_distance(data[:,0], data[:,1], 180)
        rho_q = data[:,2]

        rho = plt.figure()
        rho.add_subplot(111)
        plt.title("Density vs radius", fontsize=15)
        plt.xlabel(r"$ r \: [fm]$", fontsize=10)
        plt.ylabel(r"$ \rho_{q}$", fontsize=10) #\: {[fm}^{-3}]$
        plt.xticks(fontsize=10)
        plt.scatter(r, rho_q, 0.5) #plotting
        #plt.legend()
        plt.show()
    else:
        sys.exit('# ERROR: Cannot access file')    

def plot_current_slice(filename):
    """
    Opens the specified file, checks its validity, and creates a 2D array from
    its data. The first and second columns are x- and y-coordinates. The third,
    fourth, and fifth columns are :math:`\\vec x`, :math:`\\vec y`, :math:`\\vec z`
    current components [fm :sup:`-4`]

    Cross section distance and the total current vector :math:`\\vec j` are
    calculated. The current is then plotted against the cross section distance.

    Args:
        filename (string): name of the data set file

    Returns:
        None

    See also:
        :func:`cross_section_distance`
        :func:`current()`
    """
    filename = TXT_PATH + filename
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - x current component
        #data[:,3] - y current component
        #data[:,4] - z current component

        r = cross_section_distance(data[:,0], data[:,1], 180)
        j = vector_magnitude(data[:,2], data[:,3], data[:,4])
        # current too low - not visible on graphs

        plt.title("Current A vs radius", fontsize=15)
        plt.xlabel(r"$ r \: [fm]$", fontsize=10)
        plt.ylabel(r"$ j_{q}$", fontsize=10)
        plt.xticks(fontsize=10)
        plt.scatter(r, j, 0.5) #plotting
        #plt.legend()
        plt.show()
    else:
        sys.exit('# ERROR: Cannot access file')

def plot_pairing_field_slice(filename):
    """
    Opens the specified file, checks its validity, and creates a 2D array from
    its data. The first and second columns are x- and y-coordinates. The third
    column is the real part of the pairing field :math:`\\Delta_{rel}` [MeV]
    and the fourth is the imaginary part of the pairing field
    :math:`\\Delta_{im}` [MeV].

    Cross section distance and the total pairing field :math:`\\Delta` are
    calculated and plotted against each other.

    Args:
        filename (string): name of the data set file

    Returns:
        None

    See also:
        :func:`cross_section_distance`
        :func:`pairing_field`
    """
    filename = TXT_PATH + filename
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - delta_real
        #data[:,3] - delta_imaginary

        r = cross_section_distance(data[:,0], data[:,1], 180)
        delta, arg = pairing_field(data[:,2], data[:,3])

        plot = plt.figure()
        plot.add_subplot(111)
        plt.title("Pairing field vs radius", fontsize=15)
        plt.xlabel(r"$ r\: [fm]$", fontsize=10)
        plt.ylabel(r"$\Delta / \Delta_{\infty}$", fontsize=10)
        plt.xticks(fontsize=10)
        plt.scatter(r, delta, 0.5)
        #plt.legend()
        plt.show()
    else:
        sys.exit('# ERROR: Cannot access file')

def plot_B_q_slice(filename):
    """
    Opens the specified file, checks its validity, and creates a 2D array from
    its data. The first and second columns are x- and y-coordinates. The third
    column is the density :math:`\\rho` [fm :sup:`-3`].

    Cross section distance and the mean field potential B_q, coming from the
    variation over kinetic density, :math:`B_q`, are calculated. :math:`B_q` is
    then plotted against the cross section distance.

    Args:
        filename (string): name of the data set file

    Returns:
        None

    See also:
        :func:`cross_section_distance`
        :func:`pairing_field`
    """
    filename = TXT_PATH + filename
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - B_q_x
        #data[:,3] - B_q_y
        #data[:,4] - B_q_z

        r = cross_section_distance(data[:,0], data[:,1], 180)
        rho_n = data[:,2]
        B_q = libnest.bsk.B_q(rho_n, 0, 'n') # assuming only NeuM

        plot = plt.figure()
        plot.add_subplot(111)
        plt.title(r"$B_q$ vs radius", fontsize=15)
        plt.xlabel(r"$ r\: [fm]$", fontsize=10)
        plt.ylabel(r"$B_{q}$", fontsize=10)
        plt.xticks(fontsize=10)
        plt.scatter(r, B_q, 0.5)
        #plt.legend()
        plt.show()
    else:
        sys.exit('# ERROR: Cannot access file')

def plot_U_q_slice(filename):
    """
    Opens the specified file, checks its validity, and creates a 2D array from
    its data. The first and second columns are x- and y-coordinates. The third
    column is the real part of the pairing field :math:`\\Delta_{rel}` [MeV] and
    the fourth is the imaginary part of the pairing field :math:`\\Delta_{im}` [MeV].

    Cross section distance and the  mean field potential from density
    :math:`\\rho` variation, :math:`U_q`, are calculated. :math:`U_q` is
    plotted against the cross section distance.
    
    Args:
        filename (string): name of the data set file

    Returns:
        None

    See also:
        :func:`cross_section_distance`
        :func:`pairing_field`
    """
    filename = TXT_PATH + filename
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - U_q_x
        #data[:,3] - U_q_y
        #data[:,4] - U_q_z

        r = cross_section_distance(data[:,0], data[:,1], 180)
        rho_n = data[:,2]
        B_q = libnest.bsk.U_q(rho_n, 0, 'n') # assuming only NeuM

        plot = plt.figure()
        plot.add_subplot(111)
        plt.title(r"$U_q$ vs radius", fontsize=15)
        plt.xlabel(r"$ r\: [fm]$", fontsize=10)
        plt.ylabel(r"$ U_{q}$", fontsize=10)
        plt.xticks(fontsize=10)
        plt.scatter(r, B_q, 0.5)
        #plt.legend()
        plt.show()
    else:
        sys.exit('# ERROR: Cannot access file')

def plot_A_slice(filename):
    """
    Opens the specified file, checks its validity, and creates a 2D array from
    its data. The first and second columns are x- and y-coordinates. The third,
    fourth, and fifth columns are the x, y, and z components of the mean-field
    potential vector :math:`A`. The field is defined as the variation over thee
    components of the current j.

    Cross section distance is calculated and the resulting :math:`A` is plotted
    against it.

    Args:
        filename (string): name of the data set file

    Returns:
        None

    See also:
        :func:`cross_section_distance`
        :func:`pairing_field`
    """
    filename = TXT_PATH + filename
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - A_x
        #data[:,3] - A_y
        #data[:,4] - A_z
        #data[:,5] - ?

        r = cross_section_distance(data[:,0], data[:,1], 180)
        A = vector_magnitude(data[:,2], data[:,3], data[:,4])

        plot = plt.figure()
        plot.add_subplot(111)
        plt.title(r"$A$ vs radius", fontsize=15)
        plt.xlabel(r"$ r\: [fm]$", fontsize=10)
        plt.ylabel(r"$ A$", fontsize=10)
        plt.xticks(fontsize=10)
        plt.scatter(r, A, 0.5)
        #plt.legend()
        plt.show()
    else:
        sys.exit('# ERROR: Cannot access file')


# ================================
#           Velocities
# ================================

def plot_vsf(filename):
    filename = TXT_PATH + filename
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - rho

        r = cross_section_distance(data[:,0], data[:,1], 180)
        vsf = libnest.definitions.vsf(r)
        
        plt.figure()
        plt.title(r"$v_{sf}$ vs radius", fontsize=15)
        plt.xlabel(r"$ r\: [fm]$", fontsize=10)
        plt.ylabel(r"$ v_{sf} $ [c]", fontsize=10)
        plt.xticks(fontsize=10)
        plt.scatter(r, vsf, 0.5)
        #plt.legend()
        plt.show()
    else:
        sys.exit('# ERROR: Cannot access file')

def plot_vsf_nv(filename_density, filename_A, filename_current):
    filename_density = TXT_PATH + filename_density
    filename_A = TXT_PATH + filename_A
    filename_current = TXT_PATH + filename_current
    if file_check(filename_density) & file_check(filename_A):
        data_A = np.genfromtxt(filename_A, delimiter=' ', comments='#')
        data_A = data_A[~np.isnan(data_A).any(axis=1)]
        data_A = data_A[data_A[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - rho
        data_rho = np.genfromtxt(filename_density, delimiter=' ', comments='#')
        data_rho = data_rho[~np.isnan(data_rho).any(axis=1)]
        # data_rho = data_rho[data_rho[:, -1] != 0]
        
        data_j = np.genfromtxt(filename_current, delimiter=' ', comments='#')
        data_j = data_j[~np.isnan(data_j).any(axis=1)]
        # data_j = data_j[data_j[:, -1] != 0]

        r = cross_section_distance(data_A[:,0], data_A[:,1], 180)
        B_q = libnest.bsk.B_q(data_rho[:,2], 0., 'n')
        # print(len(data_rho[:,2]))
        # print(B_q)
        # print("B  " + str(len(B_q)))
        vsf = libnest.definitions.vsf(r)
        # print(vsf)
        # print("VSF  " + str(len(vsf)))
        A = vector_magnitude(data_A[:,2], data_A[:,3], data_A[:,4])
        # print(len(data_A[:,2]))
        # print(A)
        # print("A  " + str(len(A)))
        j = vector_magnitude(data_j[:,2], data_j[:,3], data_j[:,4])
        
        vsf_nv = libnest.definitions.vsf_NV(B_q, vsf, A)
        # print(vsf_nv)
        # print("VSF NV  " + str(len(vsf_nv)))
        v_nv = libnest.definitions.v_NV(B_q, j, data_rho[:,2], A)
        
        
        plt.figure()
        plt.title(r"$v_{sf \: NV}$ vs radius", fontsize=15)
        plt.xlabel(r"$ r\: [fm]$", fontsize=10)
        plt.ylabel(r"$ v_{sf \: NV}$ [c]", fontsize=10)
        plt.xticks(fontsize=10)
        plt.scatter(r, vsf_nv, 0.5, label="v_sf_nv")
        plt.scatter(r, v_nv, 0.5, label="v_nv")
        #plt.legend()
        plt.show()
    else:
        sys.exit('# ERROR: Cannot access file')

def plot_landau_velocity(filename_density, filename_delta):
    filename_density = TXT_PATH + filename_density
    filename_delta = TXT_PATH + filename_delta
    if file_check(filename_density) & file_check(filename_delta):
        data_delta = np.genfromtxt(filename_delta, delimiter=' ', comments='#')
        data_delta = data_delta[~np.isnan(data_delta).any(axis=1)]
        # data_delta = data_delta[data_delta[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - rho
        data_rho = np.genfromtxt(filename_density, delimiter=' ', comments='#')
        data_rho = data_rho[~np.isnan(data_rho).any(axis=1)]
        # data_rho = data_rho[data_rho[:, -1] != 0]

        r = cross_section_distance(data_rho[:,0], data_rho[:,1], 180)
        delta, arg = pairing_field(data_delta[:,2], data_delta[:,3])
        kf = libnest.definitions.rho2kf(data_rho[:,2])
        v_landau = delta / kf
        
        
        plt.figure()
        plt.title(r"$v_{sf}$ vs radius", fontsize=15)
        plt.xlabel(r"$ r\: [fm]$", fontsize=10)
        plt.ylabel(r"$ v_{sf} $ [c]", fontsize=10)
        plt.xticks(fontsize=10)
        plt.plot(r, v_landau)
        #plt.legend()
        plt.show()
    else:
        sys.exit('# ERROR: Cannot access file')

def plot_landau_critical_velocity(filename_density, filename_delta):
    filename_density = TXT_PATH + filename_density
    filename_delta = TXT_PATH + filename_delta
    if file_check(filename_density) & file_check(filename_delta):
        data_delta = np.genfromtxt(filename_delta, delimiter=' ', comments='#')
        data_delta = data_delta[~np.isnan(data_delta).any(axis=1)]
        # data_delta = data_delta[data_delta[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - rho
        data_rho = np.genfromtxt(filename_density, delimiter=' ', comments='#')
        data_rho = data_rho[~np.isnan(data_rho).any(axis=1)]
        # data_rho = data_rho[data_rho[:, -1] != 0]

        r = cross_section_distance(data_rho[:,0], data_rho[:,1], 180)
        delta, arg = pairing_field(data_delta[:,2], data_delta[:,3])
        kf = libnest.definitions.rho2kf(data_rho[:,2])
        v_landau = libnest.definitions.vLandau(delta, kf)
        v_critical = libnest.definitions.vcritical(delta, kf)
        
        
        plt.figure()
        plt.title(r"$v_{sf}$ vs radius", fontsize=15)
        plt.xlabel(r"$ r\: [fm]$", fontsize=10)
        plt.ylabel(r"$ v_{sf} $ [c]", fontsize=10)
        plt.xticks(fontsize=10)
        plt.scatter(r, v_landau, 0.5, label="Landau")
        plt.scatter(r, v_critical, 0.5, label="Critical")
        plt.legend(loc="upper right")
        plt.show()
    else:
        sys.exit('# ERROR: Cannot access file')

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
            i = np.where(np.logical_and(r>=40, r<=60))
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
            i = np.where(np.logical_and(r>=40, r<=60))
            delta, arg = pairing_field(data[i,2], data[i,3])
            delta_max_uniform.append(np.max(delta))
            temperature_uniform.append(float(file[-16:-12]))
        else:
            sys.exit('# ERROR: Cannot access file')

    # temperature_critical = delta_max/1.764
    delta_max = np.array(delta_max)
    temperature_critical = np.mean(delta_max[0:4]/1.764)
    print(temperature_critical)

    delta_theoretical_zero = 1.764*temperature_critical
    delta_theoretical_t_critical = 3.06*temperature_critical*np.sqrt(1-temperature/temperature_critical)
    #np.where((temperature>temperature_critical), 3.06*temperature_critical*np.sqrt(1-temperature/temperature_critical), delta_theoretical_t_critical)
    x = np.array(temperature)/float(temperature_critical)
    x_uniform = np.array(temperature_uniform)/float(temperature_critical)

    plt.figure()
    plt.xlim(0., 1.1)
    plt.title(r"Pairing field $\Delta$ for "+particles_nr+" particles", fontsize=15)
    plt.xlabel(r"$ T \: [MeV/k_B]$", fontsize=10)
    plt.ylabel(r"$\Delta_{max}$ [MeV]", fontsize=10)
    # plt.scatter(temperature, delta_max, label="vortex")
    # plt.scatter(temperature_uniform, delta_max_uniform, label="uniform")
    plt.plot(x, delta_max, '-o', label="vortex")
    plt.plot(x_uniform, delta_max_uniform, '-o', label="uniform")
    plt.plot(x[6:], delta_theoretical_t_critical[6:], label="delta for T~T crit", linestyle='dashed', )
    plt.legend(loc='upper right', ncol=1)
    plt.show()
    

    
if __name__ == '__main__':
    pass
