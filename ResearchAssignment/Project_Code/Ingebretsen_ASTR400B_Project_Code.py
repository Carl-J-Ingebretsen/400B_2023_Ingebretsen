#Carl Ingebretsen ASTR400B Spring 2023
#This file contains the code for my semester project in this class
#I am analyzing the merger remenant of the Milky Way and Andromeda to determine if 
#the remnant is a S0 type galaxy.

import numpy as np
from matplotlib import pyplot as plt
import astropy.units as u
import astropy.constants as const
from Readfile import Read
from CenterOfMass import CenterOfMass
from GalaxyMass import ComponentMass
from MassProfile import MassProfile

#Outline
'''
1) choose a snap shot after the merger has taken place between MW and M31.
This will be one from near the end

2)Compute the COM of the remenant

3) Plot the merger remnant stars and see if its elliptical or if it has an obvious disk component 
A disk component would make it type S0.

4) Possibly plotting a velocity dispersion would help to see if a disk is present

5) Either way then fit a Sersic profile to the galaxy assuming M/L=1
The elliptical part should be a n=4 de Vaucoulers 
If a disk is present there should be a contribution from a n~1 spiral disk like component

6) I will try to find the best fit profile for the intensity vs radius
Will one function work or maybe two?

7) I will test the code on a few different snapshots from the end of the simulation to be sure I didn't
choose an abnormal example.

'''


def main():
    '''This is the main function of this program.'''

def rotate_galaxy():
    '''Rotate the galaxy to face on'''

def fit_sersic_profile():
    '''fit a sersic profile to the galaxy'''

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
    '''Make a sersic for spirals'''