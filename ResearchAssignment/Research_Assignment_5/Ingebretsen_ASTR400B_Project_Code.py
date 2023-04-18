#Carl Ingebretsen ASTR400B Spring 2023
#This file contains the code for my semester project in this class
#I am analyzing the merger remenant of the Milky Way and Andromeda to determine if 
#the remnant is a S0 type galaxy.

import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import astropy.units as u
import astropy.constants as const
from Readfile import Read
from CenterOfMass import CenterOfMass
from GalaxyMass import ComponentMass
from MassProfile import MassProfile
from astropy.constants import G
import scipy.optimize as so

#Outline
'''
1) choose a snap shot after the merger has taken place between MW and M31.
This will be one from near the end. Probably one from 790-800

2)Compute the COM of the remenant using MW and M31 disk stars (bulge stars too?)

3) Plot the merger remnant stars and see if its elliptical or if it has an obvious disk component 
A disk component would make it type S0.

4) Possibly plotting a velocity dispersion would help to see if a disk rotation is present

5) Either way then fit a Sersic profile to the galaxy assuming M/L=1
The elliptical bulge part should be a n=4 de Vaucoulers 
If a disk is present there should be a contribution from a n~1 spiral disk like component
I will adjust until one profile or two profiles superimposed on each other fit the galaxy.

Fit with a function?

'''


def main():
    '''This is the main function of this program.'''

def rotate_galaxy():
    '''Rotate the galaxy to face on'''

def fit_sersic_profile():
    '''fit a sersic profile to the galaxy
    n=4 for the elliptical bulge
    n~1 for the spiral disk compoennt
    What to use numpy? scipy curve fit?'''

#From lab 7 this is the sersic for ellipticals
def sersicE(r,re,n,mtot):
    """This function returns the Sersic Profile for an elliptical galaxy 
    assuming the mass/light is roughly 1
    Input:
    r: (float) the distance from the center of the gal in (kpc)
    re: (float) the effective half light radius in (kpc)
    n: The (float) the sersic index
    mtot: (float) the total stellar mass of the system in Msun
    Returns:
    I: the array of floats. The surface brightness profile of elliptical Lsun/kpc^2"""
    
    #M/L is 1 so total lum is 
    lum = mtot
    #the effective surface brightness
    Ie = lum/7.2/np.pi/re**2
    #Return the surface brightness profile
    a=(r/re)**(1/n)
    b=-7.67*(a-1)
    return Ie*np.exp(b)

def seresicS(r,mtot):
    '''Make a sersic for spirals
    Set n~1 and M/L=1'''

#from lab 7
def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

#from lab 7
def density_contour(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
    """ Create a density contour plot.
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
        
    Example Usage
    -------------
     density_contour(x pos, y pos, contour res, contour res, axis, colors for contours)
     e.g.:
     density_contour(xD, yD, 80, 80, ax=ax, 
         colors=['red','orange', 'yellow', 'orange', 'yellow'])

    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))
    
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T
    fmt = {}
    
    ### Adjust Here #### 
    
    # Contour Levels Definitions
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    cont_4 = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.75))
    cont_5 = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.30))
    
    # You might need to add a few levels


    # Array of Contour levels. Adjust according to the above
    levels = [cont_5,one_sigma, cont_4,two_sigma, three_sigma][::-1]
    
    # contour level labels  Adjust accoding to the above.
    strs = ['0.30','0.68','0.75','0.98', '0.99'][::-1]

    
    ###### 
    
    if ax == None:
        contour = plt.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        plt.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)

    else:
        contour = ax.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        ax.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)
    
    return contour

#from lab 7
def RotateFrame(posI,velI):
    """a function that will rotate the position and velocity vectors
    so that the disk angular momentum is aligned with z axis. 
    
    PARAMETERS
    ----------
        posI : `array of floats`
             3D array of positions (x,y,z)
        velI : `array of floats`
             3D array of velocities (vx,vy,vz)
             
    RETURNS
    -------
        pos: `array of floats`
            rotated 3D array of positions (x,y,z) such that disk is in the XY plane
        vel: `array of floats`
            rotated 3D array of velocities (vx,vy,vz) such that disk angular momentum vector
            is in the +z direction 
    """
    
    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    # normalize the vector
    L_norm = L/np.sqrt(np.sum(L**2))


    # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
    
    # z unit vector
    z_norm = np.array([0, 0, 1])
    
    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))
    
    # dot product between L and z 
    c = np.dot(L_norm, z_norm)
    
    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

    # Rotate coordinate system
    pos = np.dot(R, posI.T).T
    vel = np.dot(R, velI.T).T
    
    return pos, vel