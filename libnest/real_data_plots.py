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
def cross_section_distance(x, y):
    """Returns the distance to the centre of the 90x90 fm box as cross-section,
    changing it from a coordinate system with origin at one corner of the box.
    
     Args:
        x (float): x-coordinate :math:`r` [fm]; 
        y (float): y=coordinate :math:`r` [fm`];

    Returns:
        float: cross section distance :math:`r` [fm]
    
    """
    r = np.sqrt((x-45)**2 + (y-45)**2)
    i = np.where(x<45)
    r[i] = -r[i]
    return r

def phi(x,y): #just in case
    """Returns the angle of the line connecting the coordinates to the centre
    of the 90x90 "box", for a coordinate system with origin at a corner of the
    box.
    
     Args:
        x (float): x-coordinate :math:`r` [fm]; 
        y (float): y=coordinate :math:`r` [fm`];

    Returns:
        float: angle :math:`\\phi` [rad]
    
    """
    return np.arctan(y-45,x-45)


# ================================
#      Real data definitions
# ================================
def pairing_field(rel, im):
    """Calculates the absolute value/modulus and the argument of the field 
    potential from its imaginary and real parts.
    
    Args:
        rel (float): real part of field potential :math:`\\Delta_{rel}` [MeV]; 
        im (float): imaginary part of field potential :math:`\\Delta_{im}` [MeV];
        
    Returns:
        float: pairing field :math:`\\Delta` [MeV];
    
    """
    return np.sqrt(im**2+rel**2), np.arctan(im,rel)

def current(x, y, z):
    """Calculates the current value from its vector x, y, and z components.
    
    Args:
        x (float): :math:`\\vec x` current component [fm :sup:`-4`]
        y (float): :math:`\\vec y` current component [fm :sup:`-4`]
        z (float): :math:`\\vec z` current component [fm :sup:`-4`]
        
    Returns:
        float: current :math:`\\vec j` [fm :sup:`-4`]
    
    """
    return np.sqrt(x**2+y**2+z**2)


# ================================
#       Plotting functions
# ================================
def plot_density(filename):
    """
    Opens the specified file, checks its validity, and creates 2D arrays from
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
        
        r = cross_section_distance(data[:,0], data[:,1])
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
    Opens the specified file, checks its validity, and creates 2D arrays from
    its data. The first and second columns are x- and y-coordinates. The third
    column is the density :math:`\\rho` [fm :sup:`-3`].
    
    A grid is created and populated using triangulation ("scipy.interpolate.triangulate")
    
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
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - delta_rel
        #data[:,3] - delta_im
        
        r = cross_section_distance(data[:,0], data[:,1])
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
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        print(data)
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - i current component
        #data[:,3] - j current component
        #data[:,4] - k current component
        
        r = cross_section_distance(data[:,0], data[:,1])
        print(r)
        j = current(data[:,2], data[:,3], data[:,4])
        print(j)

        #NOTE: treats some files as empty? to double check

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
        
        
if __name__ == '__main__':
    pass