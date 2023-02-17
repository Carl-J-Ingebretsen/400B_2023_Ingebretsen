# Homework 4: Carl Ingebretsen ASTRO400B Spring 2023
# Center of Mass Position and Velocity

# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl
from matplotlib import pyplot as plt
from astropy.constants import G

from Readfile import Read
from CenterOfMass import CenterOfMass

#define the class of massprofile
class MassProfile:
    #class description

    def __init__(self, galaxy, snap):
        '''This class constructs the mass profiles of a galaxy
        Input:
            galaxy: the name of the galaxy as a string "MW", "M31" or "M33"
            snap: the snap number as an integer
        '''
        # add a string of the filenumber to the value “000” 
        ilbl = '000' + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        #self.filename="%s "%(galaxy) + ilbl +'.txt'#This is sus
        self.filename=galaxy+'_'+ilbl+'.txt' #This may work

        #Read in the data
        self.time, self.total, self.data = Read(self.filename)
        #Assign Units to xyz?
        self.m = self.data['m']*10**10
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc

        #assign the galaxy name
        self.gname = galaxy

    def MassEnclosed(self, particle_type, radii):
        '''This function computes the mass of the galaxy enclosed inside the value of a radius.
        Inputs:
        particle_type: the integer 0,1,or 2 indicating the type of particle we want to know the mass of
        radii: an array of radius values to compute the mass at
        
        Returns:
        mass_array: an array of masses enclosed by the corresponding radii.
        '''
        #Define a center of mas
        COM = CenterOfMass(self.filename,particle_type)
        com_p = COM.COM_P(0.1)
        #initialize an array to store the mass:
        mass_array = np.zeros(len(radii))
        particle_index = np.where(self.data['type']==particle_type)[0]
        x = self.x[particle_index]
        y = self.y[particle_index]
        z = self.z[particle_index]
        #Now loop over the array summing masses
        for i in range(len(radii)):
            #Get the particles inside the radius
            r_part=((x-com_p[0])**2+(y-com_p[1])**2+(z-com_p[2])**2)**0.5
            index = np.where(r_part.value<radii[i])
            mass_inside = self.m[index]
            #Sum the masses
            mass_array[i]=np.sum(mass_inside)

        return mass_array*u.Msun


    def MassEnclosedTotal(self,radii):
        '''This function finds the total mass of the galaxy enclosed within the radius of all components
        Inputs:
        radii: an array of radii for in wich the mass should be calculated:
        Returns:
        An array of total masses for enclosed in each radii.'''
        #Buldge
        if self.gname != 'M33':
            buldge_mass=self.MassEnclosed(3,radii)
        else:
            buldge_mass=np.zeros(len(radii))
        #halo
        halo_mass=self.MassEnclosed(1,radii)
        #Disk
        disk_mass=self.MassEnclosed(2,radii)
        #Add them together but check to make sure its not empty
        return buldge_mass+halo_mass+disk_mass

    def HernquistMass(self, radius, a, Mhalo):
        '''This function calculates the halo mass for a henquist profile for a given radius and 
        scale factor a. Mhalo is the total halo mass
        Inputs:
        radius: the radius in which to calculate the mass
        a: the scale factor
        Mhalo: The total mass of the dark matter halo
        Returns:
        the mass enclosed in the Henquist profile at the radius.
        '''
        #The mass of the halo. This is the input but its here for reference.
        #i = np.where(self.data['type']==1) #Find all the halo particles
        #mass = np.sum(self.data['m'][i])*1e10*u.Msun#Total mass of the halo
        mass_array = np.zeros(len(radius))
        for i in range(len(radius)):
            mass_array[i] = Mhalo.value*radius[i].value**2/(a.value+radius[i].value)**2
        #mass = Mhalo*1e12*radius**2/(a+radius)**2*u.Msun #Hernquist mass
        return mass_array*u.Msun

    def CircularVelocity(self, particle_type, radii):
        '''This function calculates the circular speed at a radius for an enclosed mass.
        Inputs:
        particle_type: The type of particle to use for the mass
        radii: The array od radii value to calculate the speed at.
        Returns:
        an array of speeds correpsonding to each radii'''

        #Conver to new units
        G_withunits = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        #First calcualte the mass enclosed at each radius
        mass_array = self.MassEnclosed(particle_type,radii)
        #now calcualte the circualr speed
        v_array = np.zeros(len(radii))
        for i in range(len(radii)):
            v_array[i]=(G_withunits.value*mass_array[i].value/radii[i])**0.5 #Calculate the speed
        
        return np.round(v_array*u.km/u.s,2)

    def CircularVelocityTotal(self, radii):
        '''This function calculates the total circular velocity for the entire mass of the galaxy
        Input:
        radii: An array of radii for which the roation curve will be calculated.
        Output:
        v_array: The array of velocity values'''

        #Calculate the total mass of the galaxy enclosed in each radii
        #Conver to new units
        G_withunits = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        mass_tot = self.MassEnclosedTotal(radii)
        v_array = np.zeros(len(radii))
        for i in range(len(radii)):
            v_array[i]=(G_withunits.value*mass_tot[i].value/radii[i])**0.5 #Calculate the speed
        
        return np.round(v_array,2)

    def HernquistVCirc(self, radius, a, Mhalo):
        '''This function calculates the total circular velocity for the entire mass of the galaxy
        using the Henquist Profile.
        Input:
        radius: An array of radii for which the roation curve will be calculated.
        a: The scale radius in kpc
        Mhalo: The halo mass in M_sun
        Output:
        v_array: The array of velocity values'''
        #Calculate the total mass of the galaxy enclosed in each radii
        #Conver to new units
        G_withunits = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        mass = self.HernquistMass(radius,a,Mhalo)
        v_Hern=(G_withunits.value*mass.value/radius)**0.5 #Calculate the speed
        return np.round(v_Hern,2)

#test the funcitons out 
#Do this for the milky way first
'''MW = MassProfile('MW',0)
i = np.where(MW.data['type']==1) #Find all the halo particles
MW_halo_mass = np.sum(MW.data['m'][i])*1e10*u.Msun#Total mass of the halo
r = np.arange(0.1, 30.5, 1.5); print(r)
print(MW.MassEnclosed(1, r))
plt.plot(r, MW.MassEnclosed(1,r),label="Halo")#The halo
plt.plot(r, MW.MassEnclosed(2,r),label="Disk")#The disk
plt.plot(r, MW.MassEnclosed(3,r),label="Bulge")#The buldge What order?
plt.plot(r, MW.MassEnclosedTotal(r),label="Total Mass")#The total mass curve
plt.plot(r, MW.HernquistMass(r*u.kpc,1.2*u.kpc,MW_halo_mass),label="Henquist")#The henquist curve must fix
plt.title("Graph of mass enclosed for MW")
plt.xlabel("Radius in kpc")
plt.ylabel('Mass in Solar Masses')
plt.legend()
plt.semilogy()
plt.show()

#Now for M31
M31 = MassProfile('M31',0)
i = np.where(M31.data['type']==1) #Find all the halo particles
M31_halo_mass = np.sum(M31.data['m'][i])*1e10*u.Msun#Total mass of the halo
r = np.arange(0.1, 30.5, 1.5); print(r)
print(M31.MassEnclosed(1, r))
plt.plot(r, M31.MassEnclosed(1,r),label="Halo")#The halo
plt.plot(r, M31.MassEnclosed(2,r),label="Disk")#The disk
plt.plot(r, M31.MassEnclosed(3,r),label="Bulge")#The buldge What order?
plt.plot(r, M31.MassEnclosedTotal(r),label='Total Mass')#The total mass curve
plt.plot(r, M31.HernquistMass(r*u.kpc,1.5*u.kpc,M31_halo_mass),label="Henquist")#The henquist curve must fix
plt.title("Graph of mass enclosed for M31")
plt.xlabel("Radius in kpc")
plt.ylabel('Mass in Solar Masses')
plt.legend()
plt.semilogy()
plt.show()

#Now for M33: Theres a problem here
M33 = MassProfile('M33',0)
i = np.where(M33.data['type']==1) #Find all the halo particles
M33_halo_mass = np.sum(M33.data['m'][i])*1e10*u.Msun#Total mass of the halo
r = np.arange(0.1, 30.5, 1.5); print(r)
print(M33.MassEnclosed(1, r))
plt.plot(r, M33.MassEnclosed(1,r),label="Halo")#The halo
plt.plot(r, M33.MassEnclosed(2,r),label='Disk')#The disk
#plt.plot(r, M33.MassEnclosed(3,r))#The buldge What order?
plt.plot(r, M33.MassEnclosedTotal(r),label='Total Mass')#The total mass curve
plt.plot(r, M33.HernquistMass(r*u.kpc,0.75*u.kpc,M33_halo_mass),label='Henquist Profile')#The henquist curve must fix
plt.title("Graph of mass enclosed for M33")
plt.xlabel("Radius in kpc")
plt.ylabel('Mass in Solar Masses')
plt.legend()
plt.semilogy()
plt.show()

#Now for Rotation Curves
#MW = MassProfile('MW',0)
#i = np.where(MW.data['type']==1) #Find all the halo particles
#MW_halo_mass = np.sum(MW.data['m'][i])*1e10*u.Msun#Total mass of the halo
#r = np.arange(0.1, 30.5, 1.5); print(r)
#print(MW.MassEnclosed(1, r))
#r=r*u.kpc
plt.plot(r, MW.CircularVelocity(1,r),label="Halo")#The halo
plt.plot(r, MW.CircularVelocity(2,r),label="Disk")#The disk
plt.plot(r, MW.CircularVelocity(3,r),label="Bulge")#The buldge What order?
plt.plot(r, MW.CircularVelocityTotal(r),label="Total Mass")#The total mass curve
plt.plot(r, MW.HernquistVCirc(r*u.kpc,1.2*u.kpc,MW_halo_mass),label="Henquist Speed")#The henquist curve must fix
plt.title("Graph of Rotation Curve for MW")
plt.xlabel("Radius in kpc")
plt.ylabel('Velocity in km/s')
plt.legend()
plt.semilogy()
plt.show()

#Now for M31
#M31 = MassProfile('M31',0)
#i = np.where(M31.data['type']==1) #Find all the halo particles
#M31_halo_mass = np.sum(M31.data['m'][i])*1e10*u.Msun#Total mass of the halo
#r = np.arange(0.1, 30.5, 1.5); print(r)
#print(M31.MassEnclosed(1, r))
plt.plot(r, M31.CircularVelocity(1,r),label="Halo")#The halo
plt.plot(r, M31.CircularVelocity(2,r),label="Disk")#The disk
plt.plot(r, M31.CircularVelocity(3,r),label="Bulge")#The buldge What order?
plt.plot(r, M31.CircularVelocityTotal(r),label='Total Mass')#The total mass curve
plt.plot(r, M31.HernquistVCirc(r*u.kpc,1.5*u.kpc,M31_halo_mass),label="Henquist")#The henquist curve must fix
plt.title("Graph of Rotation Curve for M31")
plt.xlabel("Radius in kpc")
plt.ylabel('Velocity in km/s')
plt.legend()
plt.semilogy()
plt.show()

#Now for M33: Theres a problem here
#M33 = MassProfile('M33',0)
#i = np.where(M33.data['type']==1) #Find all the halo particles
#M33_halo_mass = np.sum(M33.data['m'][i])*1e10*u.Msun#Total mass of the halo
#r = np.arange(0.1, 30.5, 1.5); print(r)
#print(M33.MassEnclosed(1, r))
plt.plot(r, M33.CircularVelocity(1,r),label="Halo")#The halo
plt.plot(r, M33.CircularVelocity(2,r),label='Disk')#The disk
#plt.plot(r, M33.MassEnclosed(3,r))#The buldge What order?
plt.plot(r, M33.CircularVelocityTotal(r),label='Total Mass')#The total mass curve
plt.plot(r, M33.HernquistVCirc(r*u.kpc,0.75*u.kpc,M33_halo_mass),label='Henquist Profile')#The henquist curve must fix
plt.title("Graph of Rotation Curve for M33")
plt.xlabel("Radius in kpc")
plt.ylabel('Velocity in km/s')
plt.legend()
plt.semilogy()
plt.show()
'''