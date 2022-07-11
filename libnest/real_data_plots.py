# -*- coding: utf-8 -*-
"""
Plotting real data
"""
import numpy as np
import matplotlib.pyplot as plt
import libnest.definitions
import libnest.plots

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

#OPENING FILES
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

def density(filename):
    if file_check(filename):
        DATA = np.genfromtxt(filename, delimiter=' ', comments='#')
        DATA = DATA[~np.isnan(DATA).any(axis=1)]
        DATA = DATA[DATA[:, -1] != 0]
        #DATA[:,0] - x
        #DATA[:,1] - y
        #DATA[:,2] - rho_q

        #r equation!
    
        RHO = plt.figure()
        RHO.add_subplot(111)
        plt.title("Density vs radius", fontsize=15)
        plt.xlabel(r"$ r\: [fm]$", fontsize=10)
        plt.ylabel(r"$\rho \: {[fm]}^{-3}$", fontsize=10)
        plt.xticks(fontsize=10)
        plt.plot(r, DATA[:,2], linewidth=2.0, label='Fit')
        plt.legend()

        plt.show()
        
        
if __name__ == '__main__':
    pass