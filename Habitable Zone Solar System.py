# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 16:34:53 2022

@author: Aaron
"""

'''
inner solar system 
draw the goldylock zone
see if u can add another planet to the system in a stable orbit
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time 


'''
defining some constants
'''
G = 6.67e-11 
secsinyear = 3.154e7

t = np.arange(0,2 * secsinyear, secsinyear/(365/4)) #creating an array of time values from 0 to 12 years in steps of a quarter of a year

AU = 1.496e11

M0 = 1.989e30

L0 = 3.846e26

Day1 = 60*60*24

'''
creating an array to store the surface temps in
'''
Mi = 0.00001


'''
creating a class that can make any of the large objects like a planet/sun
'''

class Star:
    def __init__(self,Mass,Luminosity):
        self.x = 0
        self.y = 0
        self.Vx = 0
        self.Vy = 0
        self.mass = Mass * M0 
        self.Lum = L0 * Luminosity
        
Sun = Star(1, 1**3.5)

class SpaceObj:
    def __init__(self,y,Vx,Mass,a):
        self.x = 0
        self.y = y * AU
        self.Vx = Vx
        self.Vy = 0
        self.mass = Mass
        self.albedo = a 
        d = np.sqrt((Sun.x - self.x)**2 + (Sun.y - self.y)**2)
        self.SurfT = (Sun.Lum * (1-a)/(16*np.pi*5.67051e-8))**(1/4) * 1/np.sqrt(d) - 273
    
Mercury = SpaceObj(0.4,47870,3.3e23,0.12)
Venus = SpaceObj(0.7,35020,4.87e24,0.75)
Earth = SpaceObj(0.98,29780,5.97e24,0.3)
Mars = SpaceObj(1.5,24077,6.42e23,0.16)
    
MyObj = np.array([Sun,Mercury,Venus,Earth,Mars]) #storing all of the objects in one array

Fig1 = plt.figure(figsize = (12,7))
plt.subplots_adjust(wspace=0.6, hspace=0.8, top=0.85, left = 0.1)
plt.suptitle('Our Solar System and its habitable zone')
ax1 = Fig1.add_subplot(211)
ax2 = Fig1.add_subplot(212)


def GoldilocksZone(Star):
    iZone = np.sqrt((Star.Lum/L0)/1.1) 
    oZone = np.sqrt((Star.Lum/L0)/0.53) 
    innerC = plt.Circle((0,0),iZone, color = 'w')
    outerC = plt.Circle((0,0),oZone, color = 'c')

    ax1.add_patch(outerC)
    ax1.add_patch(innerC)
    


'''
creating arrays that will store the x,y,Vx,Vy vals of all the objects, also creating one for the masses
'''
Invals = np.array([])
Masses = np.array([])

NumofObjects = len(MyObj)
for i in range(0,len(MyObj)):
    Invals = np.append(Invals,MyObj[i].x)
    Invals = np.append(Invals,MyObj[i].y)
    Invals = np.append(Invals,MyObj[i].Vx)
    Invals = np.append(Invals,MyObj[i].Vy)
    Masses = np.append(Masses,MyObj[i].mass)

'''
defining a function to work out the acceleration
'''
def Accel(Pos1x,Pos2x,Pos1y,Pos2y,Mass):
    r1 = Pos1x - Pos2x
    r2 = Pos1y - Pos2y

    R = np.sqrt(r1**2 + r2**2)
    Acc = -G * Mass * r1 / R**3

    return Acc
    
'''
defining the function that will be used for the odeint cycle
'''
    
def Func(vals,dt):
    '''
    defining an empty array that will be filled with the new velocities and accelerations
    '''
    NewVals = np.zeros(NumofObjects * 4)
    
    '''
    B is an array of the indexs of the x pos for each object except for the ith object so that it doesnt do caluclations with itslef
    iarr is an array of the same shape as B but filled with the i value so that it is the index for the xpos of the ith object
    NewValsx and NewValsy send the x and y positions for the ith and every other object into the accel function to get a return of an array of the acceleration contributions for the ith object
    NewVals[i+2] and NewVals[i+3] add the columns of NewValsx and NewValsy to get the total accel for the ith object
    NewVals[i] and NewVals[i+1] get the x and y velocities passed into them
    '''
   
    for i in range (0,NumofObjects * 4 -1,4):
        B = np.arange(0,NumofObjects*4,4, dtype = int)
        B = np.delete(B,int(i/4))
        iarr = np.ones(shape = np.shape(B), dtype = int ) * i

        NewValsx = Accel(vals[iarr],vals[B], vals[iarr+1], vals[B+1], Masses[B//4]) 
        NewValsy = Accel(vals[iarr+1],vals[B+1], vals[iarr], vals[B], Masses[B//4])
        
        NewVals[i+2] = np.sum(NewValsx, axis = 0)
        NewVals[i+3] = np.sum(NewValsy, axis = 0)
   
        NewVals[i] = vals[i+2]
        NewVals[i+1] = vals[i+3]
           
    return NewVals

start_time = time.time()
 
'''
All the postiitons and velocities at each time step are stored in the 2D array FinalVals
'''
FinalVals = odeint(Func,Invals,t, rtol = 1e-12, atol = 1e-12)

end_time = time.time()
print('Time taken = ', end_time - start_time)

GoldilocksZone(Sun)


colours = np.array(['y*','r','m','g','b'])
colours2 = np.array(['y*','rx','mx','gx','bx'])

'''
loops in steps of 4 so that the ith value is the x pos and the i + 1 value is the y pos
i//4 is which object is being plotted so it knows what colour to plot it as
'''
for i in range(0,NumofObjects * 4 -1, 4):
    
    ax1.plot((FinalVals[:,i]- FinalVals[:,0]) / AU,(FinalVals[:,i+1] - FinalVals[:,1] )/ AU ,colours[i//4])
    if i != 0:
        ax2.plot(MyObj[i//4].mass /10e23,MyObj[i//4].SurfT,colours2[i//4])
#plt.axis('scaled')
ax1.set_xlabel('x [AU]')
ax1.set_ylabel('y [AU]')
ax2.set_xlabel('Mass [$10^{23}$ kg]')
ax2.set_ylabel('Surface Temp [deg C]')
'''
fancybox and shadow are aestetics
bbox_to_anchor moves the legend off of the plot and into the middle of the figure
'''
ax1.legend(('Sun','Mercury', 'Venus', ' Earth','Mars','Habitable Zone',''),fancybox = True, shadow = True, bbox_to_anchor=(0.95,1.25))
ax2.legend(('Mercury', 'Venus', ' Earth','Mars','Habitable Zone',''),fancybox = True, shadow = True, bbox_to_anchor=(1,1.25))



TotalE = 0
'''
calculating the change in energy, each object needs the KE but only need to work out the GPE for each pair once
'''
for i in range (0,NumofObjects * 4 -1,4):
    V = np.sqrt((FinalVals[:,i + 2])**2 + (FinalVals[:,i + 3])**2)
    TotalE += 0.5 * Masses[i//4] * V**2
    for j in range(i+4,NumofObjects * 4 - 3,4):
            R  = np.sqrt((FinalVals[:,i] - FinalVals[:,j])**2 + (FinalVals[:,i+1] - FinalVals[:,j+1])**2) 
            TotalE += -G * Masses[i//4] * Masses[j//4] / R
            
deltaE = (TotalE - TotalE[0] ) / TotalE[0]

Fig2 = plt.figure()
plt.plot(t/secsinyear, deltaE)
plt.xlabel('Time [yr]')
plt.ylabel('$\Delta E$')