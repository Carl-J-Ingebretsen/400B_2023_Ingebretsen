#Carl Ingebretsen ASTR400B Spring 2023
#This file contains the code for my semester project in this class
#I am analyzing the merger remenant of the Milky Way and Andromeda to determine if 
#the remnant is a S0 type galaxy.
#Keep disk profile to 1 and fit scale height

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

#Outline
'''
1) choose a snap shot after the merger has taken place between MW and M31.
This will be one from near the end. Snap shot 800

2)Compute the COM of the remenant using MW and M31 disk stars

3) Plot the merger remnant stars and see if its elliptical or if it has an obvious disk component 
A disk component would make it type S0.

4) Get a plot of the disk rotation velocity to check if its rotating.

5) Either way then fit a Sersic profile to the galaxy assuming M/L=1
The elliptical bulge part should be a n=4 de Vaucoulers 
If a disk is present there should be a contribution from a n~1 spiral disk like component
'''

def main():
    '''This is the main function of this program. It just calls each function that is
    being used'''

    rn,vn = get_rotated_pos_vec("Combined_800.txt") #Get rotated position and velocities
    plot_rotation('Combined',rn,vn) #Plot a velocity rotation curve
    cal_Seresic()  #Fit the light profile with a Sersic profile

def fit_sersic_profile(r_annuli,Sigma,Gal_Profile):
    '''fit a sersic profile to the galaxy
    n=4 for the elliptical bulge
    n~1 for the spiral disk component
    The function uses scipy optimize to do the curve fitting.
    
    Inputs:
    r_annuli: (array of floats) the radius values as an array
    Sigma: (array of floats) the array of luminosities of the galaxy
    Gal_profile: the Sersic profile function to be fitted

    Returns:
    params: the array of fitted parameter values
    '''

    params,pcov = curve_fit(Gal_Profile,r_annuli,Sigma) #Fit the curve
    print("Params of the fitted Seresic profile: ",params)
    print("Sersic index of the profile: ", params[2])
    #Compute 1-sigma error
    perr = np.sqrt(np.diag(pcov))
    print("The error in the Sersic index is: ", perr[2])
    return params
    

#From lab 7 this is the sersic for ellipticals
def sersicE(r,re,n,mtot):
    """This function returns the Sersic Profile for an elliptical galaxy 
    assuming the mass/light is roughly 1. The equation is from class slides
    I(r) = I_e*exp(-7.67*[(r/R_e)^(1/n)-1])

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

#Define the Sersic Profile for Spiral galaxies
def seresicS(r,re,h_r,n,mtot):
    '''Make a sersic profile for spirals for spirals
    The equation is : 
    It is taken from the slides in class slides
    Set M/L=1 under these assumptions:
    I(r) = I_0*exp(-r/h_r)+I_e*exp(-7.67*[(r/R_e)^(1/n)-1])
    
    Inputs:
    r: the array of floats for the radius
    re: the half mass/half light radius as a float
    h_r: the scale hight of the sprial disk
    n: the Sersic index
    mtot: the total mass of the galaxy
    
    Returns:
    an array of floats representing the sersic light profile'''

    #M/L is 1 so total lum is 
    lum = mtot
    #the effective surface brightness
    Ie = lum/7.2/np.pi/re**2
    I_0=Ie*10**3.33
    #Return the surface brightness profile
    a=(r/re)**(1/n)
    b=-7.67*(a-1)
    return Ie*np.exp(b)+I_0*np.exp(-r/h_r)

def seresicS_3(r,re,h_r,mtot):
    '''Make a sersic profile for spirals for spirals
    The equation is : 
    It is taken from the slides in class and the paper:
    Set M/L=1 under these assumptions
    
    Inputs:
    r: the array of floats for the radius
    re: the half mass/half light radius as a float
    h_r: the scale hight of the sprial disk
    n: the Sersic index
    mtot: the total mass of the galaxy
    
    Returns:
    an array of floats representing the sersic light profile'''
    n=1
    #M/L is 1 so total lum is 
    lum = mtot
    #the effective surface brightness
    Ie = lum/7.2/np.pi/re**2
    I_0=Ie*10**3.33
    #Return the surface brightness profile
    a=(r/re)**(1/n)
    b=-7.67*(a-1)
    return Ie*np.exp(b)+I_0*np.exp(-r/h_r)

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
    print("pdf: ", pdf)
    
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
    plt.title("Plot of the Remant velocity")
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
    plt.title("Plot of the Circular Velocity as a function of distance")
    plt.plot(R,Vcirc,color='red')
    plt.plot(-R,-Vcirc,color='red')

    # Add axis labels
    plt.xlabel('x', fontsize=22)
    plt.ylabel('vy', fontsize=22)
    plt.xlim(-100,100)

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
    '''This function will plot the galaxy face on so that it can be inspected visualy
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
    #plt.ylim(-40,40)
    #plt.xlim(-40,40)
    plt.ylim(-80,80)
    plt.xlim(-80,80)

    #adjust tick label font size
    label_size = 22
    matplotlib.rcParams['xtick.labelsize'] = label_size 
    matplotlib.rcParams['ytick.labelsize'] = label_size
    plt.show()

def cal_Seresic():
    '''Calculate the Sersic profile stuff and display it as a graph.
    This graph plots the sersic profile and fits a profile for a spiral and elliptical on top.'''
    # Create a center of mass object 
    # This lets us store the x, y, z, bulge particles 
    Comb_COM = CenterOfMass("Combined_800.txt", 2) #2 for disks 3 for bulge
    Comb_COM_p = Comb_COM.COM_P(0.1)
    x = Comb_COM.x-Comb_COM_p[0].value
    y = Comb_COM.y-Comb_COM_p[1].value
    z = Comb_COM.z-Comb_COM_p[2].value
    m = Comb_COM.m

    # Compute the surface density profile 

    # calculate the radial distances in cylindrical coordinates (and theta, too)
    cyl_r_mag = np.sqrt(x**2 + y**2) #np.sum(self.alg_r[:, :2]**2, axis=1))
    cyl_theta = np.arctan2(y,x) # self.alg_r[:, 1], self.alg_r[:, 0])

    radii = np.arange(0.1, 0.95 * cyl_r_mag.max(), 0.1)
        
    # create the mask to select particles for each radius
    # np.newaxis creates a virtual axis to make tmp_r_mag 2 dimensional
    # so that all radii can be compared simultaneously
    enc_mask = cyl_r_mag[:, np.newaxis] < np.asarray(radii).flatten()

    # calculate the enclosed masses 
    # relevant particles will be selected by enc_mask (i.e., *1)
    # outer particles will be ignored (i.e., *0)
    m_enc = np.sum(m[:, np.newaxis] * enc_mask, axis=0)

    # use the difference between nearby elements to get mass in each annulus
    m_annuli = np.diff(m_enc) # one element less then m_enc
    Sigma = m_annuli / (np.pi * (radii[1:]**2 - radii[:-1]**2))

    r_annuli = np.sqrt(radii[1:] * radii[:-1]) 

    params_spiral = fit_sersic_profile(r_annuli,Sigma,seresicS)
    params_elliptical = fit_sersic_profile(r_annuli,Sigma,sersicE)

    fig, ax = plt.subplots(figsize=(9, 8))

    # Surface Density Profile
    ax.loglog(r_annuli, Sigma, lw=2, alpha=0.8,label='Simulated Galaxy light profile') #Plot the simulated light profile
    #ax.loglog(r_annuli,seresicS(r_annuli,params[0],params[1],params[2],params[3]),lw=2, alpha=0.8,label='Fitted Galaxy with Spiral') #This works
    #ax.loglog(r_annuli,seresicS(r_annuli,params[0],params[1],params[2],m_enc),lw=2, alpha=0.8,label='Fitted Galaxy with Spiral')
    #ax.loglog(r_annuli,seresicS_3(r_annuli,params[0],params[1],params[2]),lw=2, alpha=0.8,label='Fitted Galaxy with Spiral')
    #ax.loglog(r_annuli,sersicE(r_annuli,params_E[0],params_E[1],params_E[2]),lw=2, alpha=0.8,label='Fitted Galaxy with Elliptical') #This too
    ax.loglog(r_annuli,seresicS(r_annuli,params_spiral[0],params_spiral[1],params_spiral[2],params_spiral[3]),lw=2, alpha=0.8,label='Fitted Galaxy with Spiral') #This works
    ax.loglog(r_annuli,sersicE(r_annuli,params_elliptical[0],params_elliptical[1],params_elliptical[2]),lw=2, alpha=0.8,label='Fitted Galaxy with Elliptical') #This too
    # Sersic n = 4 - de Vaucouleurs
    #plt.semilogy(r_annuli,sersicE(r_annuli,params_E[0],4,params[2]), color='red',linestyle="-.",linewidth=3, label='Sersic n=4')

    ax.set(xlabel=r"$r$ [kpc]", ylabel=r"$\Sigma$(disk par) [$10^{10} M_\odot$ / kpc$^2$]", title="Galaxy Light Profile")
    plt.savefig("Galaxy_light_profile_2")

    #set axis limits
    #plt.xlim(0.1,50)
    plt.xlim(0.1,1000)
    plt.ylim(1e-6,1)

    ax.legend(loc='best')
    fig.tight_layout()
    plt.show()

main()