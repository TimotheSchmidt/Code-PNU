import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

filename = "/Users/timotheschmidt/Documents/MINES/2A/Stage PNU/Simu/DATA_SIMS/Output00034_size0128_hyperCube.fits"
file = fits.open(filename)



def Q_tab():
    """returns the 2D array of all Q(x,y)"""
    p_0 = 1
    Bx = file[4].data
    By = file[5].data
    Bz = file[6].data
    N = file[0].data
    B  = Bx**2 + By**2 + Bz**2
    
    Tab = p_0 * N *(Bx**2 - By**2) / B 
    return np.sum(Tab, 2)



def U_tab():
    p_0 = 1
    Bx = file[4].data
    By = file[5].data
    Bz = file[6].data
    N = file[0].data
    B  = Bx**2 + By**2 + Bz**2
    
    
    Tab = 2 * p_0 * N * Bx * By /B
    return np.sum(Tab, 2)



def psi_tab():
    U, Q = U_tab(), Q_tab()
    return 0.5 * np.arctan2(U, Q)


def sigma_psi_global():
    return np.std(psi_tab())*180/np.pi

def sigma_psi_line():
    """returns an array of sigma_psi of each line"""
    return np.std(psi_tab(), 0)*180/np.pi