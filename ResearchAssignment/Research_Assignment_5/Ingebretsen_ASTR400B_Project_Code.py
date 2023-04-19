#Carl Ingebretsen ASTR400B Spring 2023
#This file contains the code for my semester project in this class
#I am analyzing the merger remenant of the Milky Way and Andromeda to determine if 
#the remnant is a S0 type galaxy.

import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import astropy.units as u
import astropy.constants as const
from ReadFile import Read
from CenterOfMass import CenterOfMass
from GalaxyMass import ComponentMass
from MassProfile import MassProfile
from astropy.constants import G
import scipy.optimize as so
import matplotlib
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit

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

    #1 Load the M31 and MW snap shots

    #2 Visualize each

    #3Try to have a combined snap shot

    #4 Visulaize it

    #5rotate combined to see rotation curve

    #6 Plot brightness profile of combined

    #M/L=1
    #7 Fit with Sereic

def fit_sersic_profile(r,Lum,L_cen):
    '''fit a sersic profile to the galaxy
    n=4 for the elliptical bulge
    n~1 for the spiral disk compoennt
    What to use numpy? scipy curve fit?'''

    #Define the combined Sersic profile as a function of two parameters for the bulge and disk
    Sers = lambda a,b,h: a*(np.e**(-(r/h)**0.25))+b*(np.e**(-(r/h)))

    #Fit with scipy optimize. Fit the a,b, and scale height

    #May need to have two scale heights
    #return Sers

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

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y)) #normed=True
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))
    
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T
    fmt = {}
    
    ### Adjust Here #### 
    '''
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
    
    return contour'''

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

def plot_rotation(name,rn,vn,snap=800):
    '''This function should plot the rotation of the galaxy
    
    Inputs:
    rn: the array of rotated positions 
    vn: the array of rotated velocities
    name: (string) name of the galaxy to be used
    snap: (int) the snap number'''

    # plot position of disk particles color coded by velocity along the 3rd axis
    # plt.scatter(pos1, pos2, c=vel1)
    # ADD HERE 
    plt.scatter(rn[:,0],rn[:,2], c=vn[:,1])

    #colorbar
    cbar = plt.colorbar()
    cbar.set_label('Vy in km/s', size=22)
    # Add axis labels
    plt.xlabel('X-direction', fontsize=22)
    plt.ylabel('Z-direction', fontsize=22)
    # calculate the 2D density of the data given
    #counts,xbins,ybins=np.histogram2d(xD[index],yD[index],bins=100,normed=LogNorm()
    #adjust tick label font size
    label_size = 22
    matplotlib.rcParams['xtick.labelsize'] = label_size 
    matplotlib.rcParams['ytick.labelsize'] = label_size

    #set axis limits
    plt.ylim(-40,40)
    plt.xlim(-40,40)
    plt.show()

    #Plot the rotation as a funtion of radius
    Gal=MassProfile(name,800) #Snap 800
    R=np.arange(0.01,40,0.1)
    Vcirc=Gal.circularVelocityTotal(R)

    # Make a phase diagram
    # MW Disk Velocity Field edge on.

    # Plot 2D Histogram one component of  Pos vs Vel 
    # ADD HERE
    plt.hist2d(rn[:,0],vn[:,1],bins=500,norm=LogNorm())
    plt.colorbar()

    # Overplot Circular Velocity from the MassProfile Code
    # ADD HERE
    plt.plot(R,Vcirc,color='red')
    plt.plot(-R,-Vcirc,color='red')

    # Add axis labels
    plt.xlabel('x', fontsize=22)
    plt.ylabel('vy', fontsize=22)

    #adjust tick label font size
    label_size = 22
    matplotlib.rcParams['xtick.labelsize'] = label_size 
    matplotlib.rcParams['ytick.labelsize'] = label_size
    plt.show()

def get_rotated_pos_vec(name,type=2):
    '''This function gets the rotated rotated frame
    Inputs:
    name: string that descrbits the name of the galaxy to be displayed
    type: 0,1,2, the type of particle to be visualized
    
    Returns:
    rn: the array of rotated positions
    vn: the array of rotated velocities'''

    # 2 for disk 1 for bulge
    COMD = CenterOfMass(name,2)
    # Compute COM of M31 using disk particles
    COMP = COMD.COM_P(0.1)
    COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])

    # Determine positions of disk particles relative to COM 
    xD = COMD.x - COMP[0].value 
    yD = COMD.y - COMP[1].value 
    zD = COMD.z - COMP[2].value 

    # total magnitude
    rtot = np.sqrt(xD**2 + yD**2 + zD**2)

    # Determine velocities of disk particles relatiev to COM motion
    vxD = COMD.vx - COMV[0].value 
    vyD = COMD.vy - COMV[1].value 
    vzD = COMD.vz - COMV[2].value 

    # total velocity 
    vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)

    # Vectors for r and v 
    r = np.array([xD,yD,zD]).T # transposed 
    v = np.array([vxD,vyD,vzD]).T

    #rotate the frame
    rn,vn = RotateFrame(r,v)
    plot_face_on_galaxy(rn)
    return rn,vn

def plot_face_on_galaxy(rn):
    '''This function will plot the galaxy face on
    Inputs:
    rn: an array of positions (floats) that describe the positions of the particles'''

    # M31 Disk Density 
    #fig, ax= plt.subplots(figsize=(10, 10))

    # plot the particle density for M31 
    # ADD HERE
    plt.hist2d(rn[:,0],rn[:,1],bins=500,norm=LogNorm(),cmap='magma')
    plt.colorbar()

    density_contour(rn[:,0],rn[:,1],80,80,colors=['magenta','green','red','yellow','gold'])

    # make the contour plot
    # x pos, y pos, contour res, contour res, axis, colors for contours.
    # ADD HERE

    # Add axis labels
    plt.xlabel('x-direction', fontsize=22)
    plt.ylabel('y-direction', fontsize=22)

    #set axis limits
    plt.ylim(-40,40)
    plt.xlim(-40,40)

    #adjust tick label font size
    label_size = 22
    matplotlib.rcParams['xtick.labelsize'] = label_size 
    matplotlib.rcParams['ytick.labelsize'] = label_size
    plt.show()

#plot_rotation_test()
#rn,vn = get_rotated_pos_vec("M31_800.txt")
#plot_rotation('M31',rn,vn)

rn,vn = get_rotated_pos_vec("Combined_800.txt")
plot_rotation('Combined',rn,vn)