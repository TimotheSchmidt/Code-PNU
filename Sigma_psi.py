import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from astropy.io import fits
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])

filename = "/Users/timotheschmidt/Documents/MINES/2A/Stage PNU/Simu/DATA_SIMS/Output00034_size0128_hyperCube.fits"
file = fits.open(filename)

    
def _Q_tab(x,y):
    p_0 = 1
    Bx = file[4].data[x,y,:]
    By = file[5].data[x,y,:]
    Bz = file[6].data[x,y,:]
    N = file[0].data[x,y,:]
    B  = Bx**2 + By**2 + Bz**2
    
    vecQ = p_0 * N *(Bx**2 - By**2) / B 
    Q = np.sum(vecQ)

    return Q

def _U_tab(x,y):
    p_0 = 1
    Bx = file[4].data[x,y,:]
    By = file[5].data[x,y,:]
    Bz = file[6].data[x,y,:]
    N = file[0].data[x,y,:]
    B  = Bx**2 + By**2 + Bz**2
    
    
    vecU = 2 * p_0 * N * Bx * By /B
    U = np.sum(vecU)

    return U

U_tab = np.vectorize(_U_tab)
Q_tab = np.vectorize(_Q_tab)

def _psi(Q,U):
    return np.arctan(U/Q)

psi = np.vectorize(_psi)

def psi_tab():
    #Q = Q_tab(np.arange(128), np.arange(128))        à essayer
    #U  = U_tab(np.arange(128), np.arange(128))
    
    Q_Tab = np.zeros((128, 128))
    U_Tab = np.zeros((128, 128))

    for x in range(128):
        for y in range(128):
             Q_Tab[x, y] = Q_tab(x, y)
             U_Tab[x, y] = U_tab(x, y)

   
    A = Q_Tab / U_Tab
    tab = 0.5 * np.arctan(A)
    return tab


def sigma_psi_global():
    return np.std(psi_tab())*180/np.pi

def sigma_psi_line():
    """returns an array of sigma_psi of each line"""
    return np.std(psi_tab(), 0)*180/np.pi