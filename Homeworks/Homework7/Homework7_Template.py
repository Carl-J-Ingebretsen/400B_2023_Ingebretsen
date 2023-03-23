
# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# Make edits where instructed - look for "****", which indicates where you need to 
# add code. 

#Carl Ingebretsen ASTR400B Spring 2023



# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
#from IPython.display import Latex


# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass import CenterOfMass
# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass

# # M33AnalyticOrbit
#now to calculate the orbit of M33



class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self,name): 
        """ This class encloses a series of functions that calcualte the orbit of M33.
        inputs:
        name: (string) the name of the file to store the output orbit in """

        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value

        ### **** store the output file name
        self.name=name
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        M33_COM = CenterOfMass("M33_000.txt",2)
        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        M33_Pos = M33_COM.COM_P(0.1,2)
        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        M33_Vel = M33_COM.COM_V(M33_Pos[0],M33_Pos[1],M33_Pos[2]).value
        M33_Pos=M33_Pos.value
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        M31_COM = CenterOfMass("M31_000.txt",2)
        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        M31_Pos = M31_COM.COM_P(0.1,2)
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        M31_Vel = M31_COM.COM_V(M31_Pos[0],M31_Pos[1],M31_Pos[2]).value
        M31_Pos=M31_Pos.value
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0 = M31_Pos-M33_Pos
        self.v0 = M31_Vel-M33_Vel
        ### get the mass of each component in M31 
        ### disk
        # **** self.rdisk = scale length (no units)
        self.rdisk = 5 #units are kpc
        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = ComponentMass("M31_000.txt", 2)*1e12
        ### bulge
        # **** self.rbulge = set scale length (no units)
        self.rbulge = 1 #units are kpc
        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = ComponentMass("M31_000.txt",3)*1e12
        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 60 #units of kpc
        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = ComponentMass("M31_000.txt",1)*1e12
    
    def HernquistAccel(self,M,r_a,r): # it is easiest if you take as an input the position VECTOR 
        """This function calculates the acceleration based on the Hernquist Profile
        
        Inputs: 
        M: (float) the Mass of the halo or bulge
        r_a: (float) the scale length
        r: (float) the relative position vector
        
        Returns:
        Hern: (floats) the acceleration vector"""
        
        ### **** Store the magnitude of the position vector
        rmag = (r[0]**2+r[1]**2+r[2]**2)**0.5
        
        ### *** Store the Acceleration
        Hern =  ((-1*self.G*M)/(rmag*(r_a+rmag)**2))*r #follow the formula in the HW instructions
        # note: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self,M,r_d,r):# it is easiest if you take as an input a position VECTOR  r 
        """ This function calculates the acceleration due to the disk of the galaxy based on the 
        formula from Miyamoto-Nagai profile 1975.
        Inputs:
        M: (float) the mass of the disk
        r_d: (float) the scale radius of the disk
        r: (float, 3 vector) the 3 vector relative position
               
        Returns:
        acc_d: (floats) the acceleration 3 vector of the disk
        """
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whole thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        R=(r[0]**2+r[1]**2)**0.5
        z_d = self.rdisk/5.0
        B=r_d+(r[2]**2-z_d**2)**0.5
        #calculate the acceleration
        z_extra = np.array([1,1,B/(r[2]**2+z_d**2)**0.5])
        acc_d = ((-1*self.G*M)/(R**2+B**2)**1.5)*r*z_extra
        return acc_d
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self,r): # input should include the position vector, r
        """This function calculates the total acceleration due to the galaxy M31
        Input:
        r: (array of floats) the relative position vector 
        
        Returns:
        acc: (array of floats) the total acceleration due to M31"""

        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        acc_bulge = self.HernquistAccel(self.Mbulge,self.rbulge,r)
        acc_halo = self.HernquistAccel(self.Mhalo,self.rhalo,r) 
        acc_disk = self.MiyamotoNagaiAccel(self.Mdisk,self.rdisk,r)   
        # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return acc_bulge+acc_halo+acc_disk
    
    
    def LeapFrog(self,dt,r,v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """This function creates a leap frog integrator. It assumes that M33 is a point mass and acceleration 
        is a pure function of r.
        Inputs:
        dt: (float) the time interval of integration
        r: (array of floats) the starting position of M33 relative to M31
        v: (array of floats) the starting velocity of M33 relative to M31
        
        Returns:
        rnew: the updated position vector
        vnew: the updated velocity vector
        """
        # predict the position at the next half timestep
        rhalf = r+v*(dt/2)
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v+self.M31Accel(rhalf)*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf+vnew*(dt/2)
        
        return rnew,vnew

    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """This function calculates the orbit of the galaxy M33
         
        Inputs:
        t0: (float) the initial time
        dt: (float) the time step
        tmax: (float) the final time to end the integration"""

        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros((int(tmax/dt)+2,7))
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t<tmax):  # as long as t has not exceeded the maximal time 
            
            # **** advance the time by one timestep, dt
            t+=dt
            # **** store the new time in the first column of the ith row
            orbit[i,0] = t
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            #r_last = np.array([orbit[i-1,1],orbit[i-1,2],orbit[i-1,3]]);v_last = np.array([orbit[i-1,4],orbit[i-1,5],orbit[i-1,6]])
            r_last = orbit[i-1,1:4];v_last=orbit[i-1,4:7]
            rnew,vnew = self.LeapFrog(dt,r_last,v_last)
    
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            orbit[i,1:4]=rnew
            orbit[i,4:7]=vnew
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i+=1
        
        
        # write the data to a file
        np.savetxt(self.name, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function

#now code the answers to the questions for this Homework
M33_orb = M33AnalyticOrbit("M33_analytical_orbit.txt") #create an instance of the object
M33_orb.OrbitIntegration(0.0,0.1,10.0) #integrate the orbit
#read in the orbit files 
M33_ana_orb = np.genfromtxt("M33_analytical_orbit.txt",names=True)
M33_hw6_orb = np.genfromtxt("Orbit_M33.txt",dtype=None,names=True)
M31_hw6_orb = np.genfromtxt("Orbit_M31.txt",dtype=None,names=True)

# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
#M31_data = np.genfromtxt("Orbit_M31.txt",dtype=None,names=True)
#print(M31_data)
#M33_data = np.genfromtxt("Orbit_M33.txt",dtype=None,names=True)
#print(M31_data)



# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  
def calculate_difference(vec_1,vec_2):
    '''This function calculates the magnitude of the differnece bewtween to given vectors
    Input:
    vec1 and vec2: two 3-vectors to find the differnece between
    Returns;
    The magnitude of the differnece between the vectors'''

    diff=np.sqrt((vec_1[0]-vec_2[0])**2+(vec_1[1]-vec_2[1])**2+(vec_1[2]-vec_2[2])**2)
    """squared_sum=0
    if len(vec_1)==len(vec_2):
        for i in range(len(vec_1)):#loop though all the terms
            squared_sum+=(vec_1[i]-vec_2[i])**2 #Calculate the squared sum
    else:
        print('invalid input')
    sum = squared_sum**0.5"""

    return diff


# Determine the magnitude of the relative position and velocities 
M31_pos = np.array([M31_hw6_orb['x'],M31_hw6_orb['y'],M31_hw6_orb['z']])
M31_vel = np.array([M31_hw6_orb['vx'],M31_hw6_orb['vy'],M31_hw6_orb['vz']])
M33_pos = np.array([M33_hw6_orb['x'],M33_hw6_orb['y'],M33_hw6_orb['z']])
M33_vel = np.array([M33_hw6_orb['vx'],M33_hw6_orb['vy'],M33_hw6_orb['vz']])
# of M33 and M31
M31_M33_sep = calculate_difference(M31_pos,M33_pos)
M31_M33_vel_sep = calculate_difference(M31_vel,M33_vel)

# Plot the Orbit of the galaxies 
#################################

plt.title("Plot of the separation of the M33 and M31")
plt.plot(M31_hw6_orb['t'],M31_M33_sep)
plt.plot(M33_ana_orb['t'],(M33_ana_orb['x']**2+M33_ana_orb['y']**2+M33_ana_orb['z']**2)**0.5)
plt.xlabel("Time (Gyrs)")
plt.ylabel("Separation (kpc)")
plt.show()

# Plot the orbital velocities of the galaxies 
#################################

plt.title("Plot of the velocity separation of the M33 and M31")
plt.plot(M31_hw6_orb['t'],M31_M33_vel_sep)
plt.plot(M33_ana_orb['t'],(M33_ana_orb['vx']**2+M33_ana_orb['vy']**2+M33_ana_orb['vz']**2)**0.5)
plt.xlabel("Time (Gyrs)")
plt.ylabel("Velocity (km/s)")
plt.show()