import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from astropy.io import fits
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])

filename = "/Users/timotheschmidt/Documents/MINES/2A/Stage PNU/Simu/DATA_SIMS/Output00034_size0128_hyperCube.fits"
file = fits.open(filename)

def Stokes_xy(file,x,y):
    n = 128
    dz = 1
    p_0 = 1
    S = [0,0]
    for z in range(n):
        rho, b_x, b_y, b_z = file[0].data[x,y,z], file[4].data[x,y,z], file[5].data[x,y,z], file[6].data[x,y,z]
        b_norm = np.sqrt(b_x**2 + b_y**2 + b_z**2)
        S[0] += rho * p_0 *(b_x**2 - b_y**2)/b_norm**2 * dz
        S[1] += 2 * rho * p_0 * b_x * b_y /b_norm**2 * dz
    return S

def psi(Q,U):
    return np.arctan(U/Q)

def psi_tab(file):
    n = 128
    tab = np.zeros((n,n))
    for x in range(n):
        for y in range(n):
            S = Stokes_xy(file, x, y)
            Q, U = S[0], S[1]
            tab[x,y] = psi(Q,U)
    return tab


def sigma_psi_global(file):
    return np.std(psi_tab(file))*180/np.pi

def sigma_psi_line(file):
    """returns an array of sigma_psi of each line"""
    return np.std(psi_tab(file), 0)*180/np.pi