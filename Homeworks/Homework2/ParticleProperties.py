#Author: Carl Ingebretsen Jan.21 2023
#Class: ASTR400B Spring 2023
#This script returns the particle properties from the galaxy simulation

#Import Read from the script readfile
from Readfile import Read
import numpy as np
import astropy.units as u

def ParticleInfo(filename, particle_type, particle_num):
    '''This function returns the properties of any selected particle from the 
    data file on the galaxy simulations.
    Inputs: filename (string), particle_type (float), particle number (int) (to identify the particle)
    Returns: Magnitude of distance in  kpc (round to 3 decimal places)
    Magnitude of distance in  lyrs (round to 3 decimal places), Magnitude of velocity in km\s
    (round to 3 decimal places), mass in units of M_sol. '''

    #Read in the file again
    time, tot_part_num, data = Read(filename)
    #Find the particles of the type asked for
    index = np.where(data['type']==particle_type)#Get the indices of all of the same type.
     #Split all the data up into its respective components:
    x = data['x'][index]*u.kpc
    y = data['y'][index]*u.kpc
    z = data['z'][index]*u.kpc
    #Speeds in each direction
    vx = data['vx'][index]*u.kilometer/u.second
    vy = data['vy'][index]*u.kilometer/u.second
    vz = data['vz'][index]*u.kilometer/u.second
    #The masses
    m = (data['m'][index])*1e10
    #Calculate the quantities for the particle of interest.
    #Mass in solar masses
    mass=m[particle_num]*u.Msun
    #Magnitude of distance
    d = (x[particle_num]**2+y[particle_num]**2+z[particle_num]**2)**0.5
    #in light years
    d_lyrs = np.around(d.to(u.lyr),3)
    #in kpc at 3 decimal places
    d = np.around(d,3)
    #Magnitude of velocity
    v = np.around((vx[particle_num]**2+vy[particle_num]**2+vz[particle_num]**2)**0.5,3)

    #test the function
    '''print("mass: ", mass)
    print("distance (lyrs): ", d_lyrs)
    print('distance (kpc): ', d)
    print("Velocity: ", v)'''
    #Return all the results
    return d,d_lyrs,v,mass

#test the function
#ParticleInfo("MW_000.txt", 2.0, 99)

