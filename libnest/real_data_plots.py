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
    Raises: FileNotFoundError
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
    changing it from a coordinate system with origin at one end of the box"""
    r = np.sqrt((x-45)**2 + (y-45)**2)
    i = np.where(x<45)
    r[i] = -r[i]
    return r

def phi(x,y): #just in case
    """Returns phi for the polar coordinates"""
    return np.arctan(y-45,x-45)


# ================================
#      Real data definitions
# ================================
def pairing_field(rel, im):
    """Returns the absolute value/modulus and the argument of the field 
    potential from its imaginary and real parts"""
    return np.sqrt(im**2+rel**2), np.arctan(im,rel)

def current(i, j, k):
    """Returns the current value from its i, j, and k components"""
    return np.sqrt(i**2+j**2+k**2)


# ================================
#       Plotting functions
# ================================
def plot_density(filename):
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
    if file_check(filename):
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        data = data[~np.isnan(data).any(axis=1)]
        data = data[data[:, -1] != 0]
        #data[:,0] - x
        #data[:,1] - y
        #data[:,2] - rho_q
        # x in one direction, y in another??
        
        x = data[:,0] - 45
        y = data[:,1] - 45
        
        xx, yy = np.meshgrid(x, y)
        rho_bulk = 0.00590448
        rho_ratio = data[:,2]/rho_bulk
        
        
        #interpolation by triangulation
        zi = scipy.interpolate.griddata((x, y), rho_ratio, (xx, yy), method='nearest')
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