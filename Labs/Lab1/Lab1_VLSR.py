
# # In Class Lab 1
# Must be uploaded to your Github repository under a "Labs/Lab1" folder by 5 PM on Jan 31st 2023

# ## Part A:  The Local Standard of Rest
# Proper motion of Sgr A* from Reid & Brunthaler 2004
# $\mu = 6.379$ mas/yr 
# 
# Peculiar motion of the sun, $v_\odot$ = 12.24 km/s  (Schonrich 2010)
# 
# 
# $v_{tan} = 4.74 \frac{\mu}{\rm mas/yr} \frac{R_o}{\rm kpc} = V_{LSR} + v_\odot$
# 
# 
# ### a)
# 
# Create a function called VLSR to compute the local standard of res (V$_{LSR}$).
# 
# The function should take as input: the solar radius (R$_o$), the proper motion (mu)
# and the peculiar motion of the sun in the $v_\odot$ direction.
# 
# Compute V$_{LSR}$ using three different values R$_o$: 
# 1. Water Maser Distance for the Sun :  R$_o$ = 8.34 kpc   (Reid 2014 ApJ 783) 
# 2. GRAVITY Collaboration Distance for the Sun:  R$_o$ = 8.178 kpc   (Abuter+2019 A&A 625)
# 3. Value for Distance to Sun listed in Sparke & Gallagher : R$_o$ = 7.9 kpc 
# 



# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from astropy import constants as const # import astropy constants


def VLSR(Sol_rad, mu=6.379, pec_motion=12.24*u.km/u.s):
    '''This function calculates the local standard of rest V_LSR.
    The inputs are:
    the radius the sun is at in the galaxy (kpc), 
    the proper motion mu in (mas/yr) mu=6.379 mas/yr (Reid & Brunthaler 2004)
    and the peculiar motion in km/s 12.24 km/s (SChonrich 2010).
    
    It returns the local standard of rest velocity in km/s'''
    return 4.74*mu*Sol_rad/u.kpc*u.km/u.s-pec_motion

#compute the motion
#V peculiar is 12.24 km/s
#Using water maser R_sol 8.34 kpc
v_pec = 12.24*u.km/u.s
mu = 6.379#*u.uas #mas?? milli or micro
R_Reid = 8.34*u.kpc #Reid 2014
VLSR_Reid=VLSR(R_Reid)
print("VLSR computed with water maser: ", VLSR(R_Reid))
R_Gravity=8.178*u.kpc #Abuter+2019
VLSR_Gravity = VLSR(R_Gravity)
print("VLSR computed with gravity collab value: ", VLSR(R_Gravity))
R_Sparke=7.9*u.kpc
VLSR_Sparke = VLSR(R_Sparke)
print("VLSR computed with Sparke and galleger value: ", VLSR(R_Sparke))
#VLSR computed with water maser:  239.9320764 km / s
#VLSR computed with gravity collab value:  235.03376988000002 km / s
#VLSR computed with Sparke and galleger value:  226.628034 km / s

# ### b)
# 
# compute the orbital period of the sun in Gyr using R$_o$ from the GRAVITY Collaboration (assume circular orbit)
# 
# Note that 1 km/s $\sim$ 1kpc/Gyr
def TorbSun(R,v):
    '''This function calculates the orbital period of the sun about the galaxy.
    T = 2*pi*R_sun /v
    
    Inputs: R: The radius of the sun orbit in kpc (astropy)
            v: The astropy quantity velocity of the sun in km/s
            
    Output: The orbital period of the sun in Gyr'''
    VkpcGyr = v.to(u.kpc/u.Gyr) #conver v from km/s it kpc/Gyr
    return 2*np.pi*R/VkpcGyr

#Use gravity collab result for VLSR
VSun = VLSR_Gravity+v_pec
T_Sun_orb =  TorbSun(R_Gravity,VSun)
print("Orbital period of the sun: ", T_Sun_orb)
#Orbital period of the sun:  0.20318680562272234 Gyr

# ### c)
# 
# Compute the number of rotations about the GC over the age of the universe (13.8 Gyr)
#Age 
Universe_age = 13.8*u.Gyr
print("Number of times the sun orbited the galaxy: ", Universe_age/T_Sun_orb)



# ## Part B  Dark Matter Density Profiles
# 
# ### a)
# Try out Fitting Rotation Curves 
# [here](http://wittman.physics.ucdavis.edu/Animations/RotationCurve/GalacticRotation.html)
# Point Mass =0
# Star mass 78*10**9 Msol
# Dark Matter mass 12*10**11 Msol
# 
# ### b)
# 
# 
# In the Isothermal Sphere model, what is the mass enclosed within the solar radius (R$_o$) in units of M$_\odot$? 
# 
# Recall that for the Isothermal sphere :
# $\rho(r) = \frac{V_{LSR}^2}{4\pi G r^2}$
# 
# Where $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$, r is in kpc and $V_{LSR}$ is in km/s
# 
# What about at 260 kpc (in units of  M$_\odot$) ? 

#density is rho = VLSR^2 / (4*pi*G*R^2)
#need to integral rho over volume
#rho*4*pi*r^2*dr=VLSR**2/(4*pi*G*r**2) *4*pi*r**2*dr
#=VLSR**2/G *dr

#gravitational constant
Grav = const.G.to(u.kpc**3/u.Gyr**2/u.Msun)

def MassIso(r, VLSR):
    '''this function will compute the dark matter mass enclosed within a given
    distnace asuming an iso thermal sphere model
    M=VLSR**2 / G * r

    Inputs:
    r: (astropy unit kpc) distance to galactic center
    VLSR: (astropy km/s) the velocity of local standard of rest
    
    Returns: the enclosed dark matter mass in M_sun (astropy)'''
    #convert km/s to kpc/Gyr
    VkpcGyr_2 = VLSR.to(u.kpc/u.Gyr) #conver v from km/s it kpc/Gyr
    M=VkpcGyr_2**2/Grav*r
    return M

M_enclosed_solar=MassIso(R_Gravity, VLSR_Gravity)
print("The mass enclosed within the solar radius: ", M_enclosed_solar)
print(f"{M_enclosed_solar:.2e}")

#Compute again for r=260kpc
M_enclosed_260=MassIso(260*u.kpc, VLSR_Gravity)
print("The mass enclosed within 260 kpc: ")
print(f"{M_enclosed_260:.2e}")

# ## c) 
# 
# The Leo I satellite is one of the fastest moving satellite galaxies we know. 
# 
# 
# It is moving with 3D velocity of magnitude: Vtot = 196 km/s at a distance of 260 kpc (Sohn 2013 ApJ 768)
# 
# If we assume that Leo I is moving at the escape speed:
# 
# $v_{esc}^2 = 2|\Phi| = 2 \int G \frac{\rho(r)}{r}dV $ 
# 
# and assuming the Milky Way is well modeled by a Hernquist Sphere with a scale radius of $a$= 30 kpc, what is the minimum mass of the Milky Way (in units of M$_\odot$) ?  
# 
# How does this compare to estimates of the mass assuming the Isothermal Sphere model at 260 kpc (from your answer above)

#The speed of the Leo I dwarf galaxy and its distance are:
v_Leo = 196*u.km/u.s
r_Leo = 260*u.kpc #Both Sohn 2013

def Hernquist_Sphere(v,r,a=30*u.kpc):
    """This function calculates the enclosed mass of the galaxy at radius r
    assuming a Hernquist Sphere model with a scale radius of a.
    The Hernquist modle is rho(r)=M/(2pi)*a/(r(r+a)^3)
    The velocity of the galaxy is:
    v^2=2 int_0^r[Grho(r)/r dV]. This results in an expression for the mass of:
    M=v^2a/(2G)*(r^2+2ar+a^2)/(r(r+2a)).
    
    Inputs:
    v: the speed of the galaxy in km/s (astropy)
    r: the radius in which the mass is to be determined (kpc)
    a: the scale radius of the model in (kpc)
    Returns:
    M: the enclosed mass in M_sol."""
    
    v_kpcGyr = v.to(u.kpc/u.Gyr) #conver v from km/s it kpc/Gyr
    M=(v_kpcGyr**2*a)/(2*Grav)*(r**2+2*r*a+a**2)/(r*(r+a))
    return M

#Mass determined with hernquist sphere:
M_Hern = Hernquist_Sphere(v_Leo,r_Leo)
print("The mass with the Hernquist Sphere model is: ", M_Hern)
print(f"{M_Hern:.2e}")
#mass with isothermal spehre:
M_iso_Leo = MassIso(r_Leo,v_Leo)
print("The mass with the isothermal Sphere model is: ", M_iso_Leo)
print(f"{M_iso_Leo:.2e}")
print("The Hernquist Model gives a mass that is about a factor of 10 lower than\
      the isothermal sphere model for the same radius and speed.")



