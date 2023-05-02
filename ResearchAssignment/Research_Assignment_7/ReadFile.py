#Author: Carl Ingebretsen 
#Class: ASTR400B Spring 2023
#This is the ReadFile program for HW2. This program will open and read the simulation files given in this class.

import numpy as np
import astropy.units as u

def Read(filename):
    '''This is the read file function that reads in the galaxy simualtion file.
        Input: the filename as a string
        Output: the time as a float (Myr), the tota number of particles as a integer, and the particle type (1=Dark matter, 
        2=disk stars, 3=halo stars), mass (10^10 M_sol)
        x,y,z,(kpc), vx, vy, vz (km/s) as an array. '''

    file = open(filename) #Open the file that was passed into the function
    line_1 = file.readline()#Read the first line (which is the time in Myr)
    label, val = line_1.split()#Split the line at the space into its description and value
    time = float(val)*u.Myr #Conver the time to a float and add units of Myr

    line_2 = file.readline()#Read the first line (which is the total particles)
    label_2, val_2 = line_2.split()#Split the line at the space into its description and value
    tot_particles = int(val_2) #Conver the total numebr of particles to a integer
    file.close()#Close the file

    #Store the rest of the file as an array
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)

    #return all the results from the file
    return time, tot_particles, data

#Test the function
#I kept these commands but commented them out in case I need to retest it in the future.
'''file_name="MW_000.txt"
t,tot,d = Read(file_name)
print("time: ", t)
print("total number of particles: ", tot)
print("data: ", d['x'][1])'''
