# -*- coding: utf-8 -*-
"""

Created on Fri Apr  3 14:46:48 2020

@author: viktor djurberg
Modules containg varibles about the system and the material. 

"""
import numpy as np

# sample information 
length=300e-6      #sample thickness (or large value for "infinite" thickness)
width= 300e-6       # sample width      (or large value for "infinite" widh)
incycles =100    #number of electrons (should be multiple of 6 for equal valley distribution)
bins_z = 100     #number of bis to divide Lenght
bins_xy = 29     #number of bis to divide width
bins_t=100    #number of temporal bins to record current data in 
min_t=-2e-9  # first travel time
max_t= 180e-9 # maximum travel time # ?set to stop_t ?
T_L = 78    # lattice temprature K
B = 0.0     # magnetic field 
U = 15      # Voltage 


#Natural and material constants
#Constans for diamond 
me=9.10953e-31 #electron mass in kg
kB=1.3807e-23  #Bolzmann constant J/K
q=1.602e-19    #elementary charge As
ml=1.56*me     #effective masses, values from N. Naka cyclotron paper. 
mt=0.28*me
m = np.array(( (ml, mt, mt), (ml, mt, mt), (mt, ml, mt), (mt, ml, mt), (mt, mt, ml), (mt, mt, ml) ),float) #array with masses in (x,y,z) direction for each valley
Xid = -6.57*q  #Dilatation deformation potential (Joule) #this value gives good fit
Xiu = 16.5*q   #Uniaxial deformation potential (Joule) #this value gives good fit  


#Advanced options 
intervalley_scattering = 0 #0 = intervalley scattering OFF, 1= intervalley scattering ON
auto_gamma_correct = 3000 #1000 #3000 # # scatterings to correct value of gamma, =0 means autocorrect off
max_scats=15000000     #maximum number of scattering events 
stop_t=1000.0e-9       #max time for simulation in seconds
              
