#Author: Carl Ingebretsen, 01/29/2023 ASTR400B
#This program contains a function that will return the total mass
#of a compnent of given galaxy. the components being either the 
#disk mass, buldge mass or dark matter mass.

from Readfile import Read
import numpy as np
import astropy.units as u

def ComponentMass(filename, particle_type):
    '''This function sums the total mass of a compoent of the galaxy.
    The componets are either (1) Halo (2) Disk or (3) Buldge.
    
    Inputs:
    filename: The name of the file as a string that the galaxy is saved in
    particle_type: an integer of value 1,2 or 3 denoting the component to be summed.
    
    Output:
    The total mass of that component of the galaxy in units of 
    10^12 solar masses.'''

    time_st, tot_particles, data = Read(filename) #Read the file and return the data
    #Now sum the mass of the component
    i = np.where(data['type']==particle_type) #Find all the indices of the desired type
    mass = np.sum(data['m'][i])*1e10*u.Msun#Sum the mass (masses where in 10^10 sol mass)
    mass = np.round(mass,3)
    return mass

#Test it out: Compute the mass of each component, the total mass and the baryon fraction
'''m_31_disk = ComponentMass('M31_000.txt',2)
print("M31 disk mass is: ", f"{m_31_disk:.3e}")
m_31_halo = ComponentMass('M31_000.txt',1) #Halo is the dark matter
print("M31 halo mass is: ", f"{m_31_halo:.3e}")
m_31_buldge = ComponentMass('M31_000.txt',3)
print("M31 buldge mass is: ", f"{m_31_buldge:.3e}")
m_31_tot = m_31_disk+m_31_buldge+m_31_halo #total mass
print("The total mass of M31 is: ", f"{m_31_tot:.3e}")
f_m_31 = (m_31_buldge+m_31_disk)/m_31_tot #baryon fraction
print("The baryon fraction of M_31 is: ", f"{f_m_31:.3e}")

m_33_disk = ComponentMass('M33_000.txt',2)
print("M33 disk mass is: ", f"{m_33_disk:.3e}")
m_33_halo = ComponentMass('M33_000.txt',1)
print("M33 halo mass is: ", f"{m_33_halo:.3e}")
m_33_tot = m_33_disk+m_33_halo#total mass
print("The total mass of M33 is: ", f"{m_33_tot:.3e}")
f_m_33 = (m_33_disk)/m_33_tot #baryon fraction
print("The baryon fraction of M_33 is: ", f"{f_m_33:.3e}")

MW_disk = ComponentMass('MW_000.txt',2)
print("MW disk mass is: ", f"{MW_disk:.3e}")
MW_halo = ComponentMass('MW_000.txt',1)
print("MW halo mass is: ", f"{MW_halo:.3e}")
MW_buldge = ComponentMass('MW_000.txt',3)
print("MW buldge mass is: ", f"{MW_buldge:.3e}")
MW_tot = MW_disk+MW_buldge+MW_halo#total mass
print("The total mass of MW is: ", f"{MW_tot:.3e}")
f_MW = (MW_buldge+MW_disk)/MW_tot #baryon fraction
print("The baryon fraction of MW is: ", f"{f_MW:.3e}")

#For the whole local group
local_group_mass = m_31_tot+m_33_tot+MW_tot
print("The total mass of the local group is: ", f"{local_group_mass:.3e}")
f_local_group = (m_31_buldge+m_31_disk+m_33_disk+MW_buldge+MW_disk)/local_group_mass
print("The baryon fraction of the local group is: ", f"{f_local_group:.3e}")'''

