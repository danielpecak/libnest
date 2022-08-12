# -*- coding: utf-8 -*-
"""
Plotting real data
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import libnest.definitions
import libnest.plots
import scipy.interpolate

#FILENAMES:
#N135i.2_current.txt
#N135i.2_density.txt
#N960i.3_current.txt
#N960i.3_density.txt
#N2600i_A.txt
#N2600i_current.txt
#N2600i_delta.txt
#N2600i_density.txt
#N2600i_info.txt
#N4160.txt
#N5000i_current.txt
#N5000i_density.txt
#N8500i_current.txt
#N8500i_density.txt
#N15000i_current.txt
#N15000i_density.txt

TXT_PATH = "C:\\Users\\aleks\\OneDrive\\Dokumenty\\libnest\\txt\\"
#depends on the user (path did not work in main)

# ================================
#         Opening Files
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


# ================================
#           Coordinates
# ================================
def cross_section_distance(x, y, size):
    """
    Returns the distance to the centre of the 90x90 fm box as cross-section,
    changing it from a coordinate system with origin at one corner of the box.
    
     Args:
        x (float): x-coordinate of :math:`r` [fm]; 
        y (float): y=coordinate of :math:`r` [fm`];
        size (float): full size of one axis of the square/box [fm`]

    Returns:
        float: cross section distance :math:`r` [fm]
    
    """
    n = size/2
    r = np.sqrt((x-n)**2 + (y-n)**2)
    i = np.where(x<n)
    r[i] = -r[i]
    return r

def phi(x,y, size): #just in case
    """
    Returns the angle of the line connecting the coordinates to the centre
    of the 90x90 "box", for a coordinate system with origin at a corner of the
    box.
    
     Args:
        x (float): x-coordinate :math:`r` [fm]; 
        y (float): y=coordinate :math:`r` [fm`];
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
        rel (float): real part of field potential :math:`\\Delta_{rel}` [MeV]; 
        im (float): imaginary part of field potential :math:`\\Delta_{im}` [MeV];
        
    Returns:
        float: pairing field :math:`\\Delta` [MeV];
    
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
        #print(data[:,0]-45)
        
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
#size is actually 180x180 for all, not 120x120
#possibly not 24 slices but 81 (see end of the function below)
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

    rho_total = np.sum(rho_q)
    print(rho_total) #gives 2.66 instead of 9 or 216 for N216_T0.00.2_density.txt etc
    print(rho_total*24) #generally fits with *81 instead of *24 slices for all
    print(216/rho_total)
    
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
# does not seem to form a clear pattern
# what is the last column in the txt files? not described in legend
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


if __name__ == '__main__':
    pass