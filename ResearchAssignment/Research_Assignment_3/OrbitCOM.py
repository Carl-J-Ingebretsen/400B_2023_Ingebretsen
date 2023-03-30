#Carl Ingebretsen Howork 6 Feb 23, 2023

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass import CenterOfMass




def OrbitCOM(galaxy,start,end,n):
    '''This function computes the time  and COM position and velocity
    vectors of a galaxy for in each snapshot and saves to a file.
    
    Inputs:
        galaxy: the name of the galaxy "MW" etc.
        start: the number of the first snapshot to be read in
        end: the number of the last snapshot
        n: the integer indicating the intervals over which you will return COM
    '''
    
    # compose the filename for output
    fileout = "Orbit_"+galaxy+".txt"
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1 #For the COM class
    volDec=2
    if galaxy=="M31":
        volDec = 4
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    if n!=0 or (start<end):#(end<start)
        snap_ids = np.arange(start,end,step=n)#Actually n as the step size?
    else:
        print("Invalid input")
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids),7]) #check to be sure
    
    # a for loop 
    for  i, snap_id in enumerate(snap_ids):# loop over files
        
        # compose the data filename (be careful about the folder)
        ilbl = '000' + str(snap_id)#make the file name to read in
        ilbl = ilbl[-3:]
        filename=galaxy+'_'+ilbl+'.txt'#This is the filename
        # Initialize an instance of CenterOfMass class, using disk particles
        #Create a COM object for the disk particles
        disk_COM = CenterOfMass(filename, 2)#Check number for disk 1 or 2
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        com_pos = disk_COM.COM_P(delta,volDec) #position COM
        com_vel = disk_COM.COM_V(com_pos[0],com_pos[1],com_pos[2])#vel COM
    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value)
        orbit[i,0]=disk_COM.time.value/1000
        orbit[i,1]=com_pos[0].value
        orbit[i,2]=com_pos[1].value
        orbit[i,3]=com_pos[2].value
        orbit[i,4]=com_vel[0].value
        orbit[i,5]=com_vel[1].value
        orbit[i,6]=com_vel[2].value
        # note that you can store 
        # a[i] = var1, *tuple(array1)

        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))




# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 

#Test first
OrbitCOM("MW",0,800,5)
OrbitCOM("M31",0,800,5)
OrbitCOM("M33",0,800,5)
#Apparently it works so far
#Generate the files
#OrbitCOM("MW",0,800,5)#Number of shots is 160? n= 5 doesn't do what I want
#OrbitCOM("M31",0,800,5)#Or maybe it does
#OrbitCOM("M33",0,800,5)

# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
MW_data = np.genfromtxt("Orbit_MW.txt",dtype=None,names=True)
print(MW_data)
M31_data = np.genfromtxt("Orbit_M31.txt",dtype=None,names=True)
print(M31_data)
M33_data = np.genfromtxt("Orbit_M33.txt",dtype=None,names=True)
print(M31_data)



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
MW_pos = np.array([MW_data['x'],MW_data['y'],MW_data['z']])
MW_vel = np.array([MW_data['vx'],MW_data['vy'],MW_data['vz']])
M31_pos = np.array([M31_data['x'],M31_data['y'],M31_data['z']])
M31_vel = np.array([M31_data['vx'],M31_data['vy'],M31_data['vz']])
M33_pos = np.array([M33_data['x'],M33_data['y'],M33_data['z']])
M33_vel = np.array([M33_data['vx'],M33_data['vy'],M33_data['vz']])
# of MW and M31
MW_M31_sep = calculate_difference(MW_pos,M31_pos)
MW_M31_vel_sep = calculate_difference(MW_vel,M31_vel)
# of M33 and M31
M31_M33_sep = calculate_difference(M31_pos,M33_pos)
M31_M33_vel_sep = calculate_difference(M31_vel,M33_vel)

# Plot the Orbit of the galaxies 
#################################
plt.title("Plot of the separation of the Milky Way and M31")
plt.plot(MW_data['t'],MW_M31_sep)
plt.xlabel("Time (Gyrs)")
plt.ylabel("Separation (kpc)")
plt.show()

plt.title("Plot of the separation of the M33 and M31")
plt.plot(M31_data['t'],M31_M33_sep)
plt.xlabel("Time (Gyrs)")
plt.ylabel("Separation (kpc)")
plt.show()

# Plot the orbital velocities of the galaxies 
#################################

plt.title("Plot of the separation of the Milky Way and M31")
plt.plot(MW_data['t'],MW_M31_vel_sep)
plt.xlabel("Time (Gyrs)")
plt.ylabel("Velocity (km/s)")
plt.show()

plt.title("Plot of the velocity separation of the M33 and M31")
plt.plot(M31_data['t'],M31_M33_vel_sep)
plt.xlabel("Time (Gyrs)")
plt.ylabel("Velocity (km/s)")
plt.show()