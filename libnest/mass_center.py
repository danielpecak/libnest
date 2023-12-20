import numpy as np
import matplotlib.pyplot as plt

def plot_mass_center(data):
    [nx,ny, nz] = [data.xyz[i].size for i in range(3)]
    pos_x=[]
    pos_y=[]
    pos_z=[]
    time=[]
    for t in range(data.rho_n.shape[0]):
        mass=data.rho_p[t]+data.rho_n[t]-data.rho_n[t][0, 0, 0]
        total_mass=np.sum(mass)
        x=0
        y=0
        z=0
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x+=data.xyz[0][i,0,0]*mass[i, j, k]
                    y+=data.xyz[1][0,j,0]*mass[i, j, k]
                    z+=data.xyz[2][0,0,k]*mass[i, j, k]
        time.append(t)
        pos_x.append(x/total_mass)
        pos_y.append(y/total_mass)
        pos_z.append(z/total_mass)
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(10, 16))    
    ax = ax.flatten()
    i=0
    ax[i].plot(time, pos_x,'-',label="Mass Center on x axis")
    i+=1
    ax[i].plot(time, pos_y,'-',label="Mass Center on y axis")
    i+=1
    ax[i].plot(time, pos_z,'-',label="Mass Center on z axis")
    labels=["x", "y", "z"]
    for i in range(3):
        ax[i].set_ylabel(labels[i])
        ax[i].legend(loc='upper left')
        ax[i].set_xlabel("time")
    plt.savefig("mass_center.png")
    plt.close()