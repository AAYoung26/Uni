# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 13:02:45 2022

@author: Aaron
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time 


fig = plt.figure(figsize = (12,8))
plt.subplots_adjust(wspace=0.6, hspace=0.8, top=0.85, left = 0.15)
plt.suptitle('Effect on the outer solar system when a large body passes by')

Titles = np.array(['New Star is 2 Solar Masses','New Star is 4 Solar Masses','New Star is 6 Solar Masses','New Star is 8 Solar Masses','New Star is 10 Solar Masses','New Star is 12 Solar Masses','New Star is 14 Solar Masses','New Star is 16 Solar Masses'])

'''
defining some constants
'''
G = 6.67 * 10**(-11) 
secsinyear = 3.154*10**7

t = np.arange(0,100 * secsinyear, secsinyear/(365/4)) #creating an array of time values from 0 to 100 years in steps of a quarter of a year

AU = 1.496 * 10**11

M0 = 1.989*10**30

Day1 = 60*60*24


'''
creating a class that can make any of the large objects like a planet/sun
'''
class SpaceObj:
    def __init__(self,x,y,Period,Mass):
        self.x = x * AU
        self.y = y * AU
        self.Vx = -225000
        self.Vy = self.x * 2 * np.pi / (Period * Day1)
        self.mass = Mass
    

Sun = SpaceObj(0,0,1,1 * M0)
Jupiter = SpaceObj(5.2,0,4331,1.9e27)
Saturn = SpaceObj(9.5,0,10747,5.68e26)
Uranus = SpaceObj(19,0,30589,8.68e25)
Neptune = SpaceObj(30,0,59800,1.02e26)
Pluto = SpaceObj(39.5,0,90560,1.31e22)
Star = SpaceObj(100,69.5,1,0.5*M0)
Star.Vy = 0
Star.Vx -= 10000



MyObj = np.array([Sun,Jupiter,Saturn,Uranus,Neptune,Pluto,Star]) #storing all of the objects in one array

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
starting a loop that will increase the mass of the asteroid each time
'''
for ma in range (1,9):
    Masses[6] = ma * 2 * M0
    '''
    defining a function to work out the acceleration
    '''

    def Accel(Pos1x,Pos2x,Pos1y,Pos2y,Mass):
        r1 = Pos1x - Pos2x
        r2 = Pos1y - Pos2y
    
        R = np.sqrt(r1**2 + r2**2)
        if R == 0 :
            R = 0.001
        Acc = -G * Mass * r1 / R**3
    
        return Acc
    
    '''
    defining the function that will be used for the odeint cycle
    the first block of code is a vectorised verison of the bottom block of code but it is slower for this number of bodies
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
  
       
        # for i in range (0,NumofObjects * 4 -1,4):
        #     B = np.arange(0,NumofObjects*4,4, dtype = int)
        #     B = np.delete(B,int(i/4))
        #     iarr = np.ones(shape = np.shape(B), dtype = int ) * i
    
        #     NewValsx = Accel(vals[iarr],vals[B], vals[iarr+1], vals[B+1], Masses[B//4])   #vectorised but it is slower for this number of bodies
        #     NewValsy = Accel(vals[iarr+1],vals[B+1], vals[iarr], vals[B], Masses[B//4])
            
        #     NewVals[i+2] = np.sum(NewValsx, axis = 0)
        #     NewVals[i+3] = np.sum(NewValsy, axis = 0)
       
        #     NewVals[i] = vals[i+2]
        #     NewVals[i+1] = vals[i+3]
            
            
        '''
        the outer loop gets the first object, the inner loop then works out the acceleration contribution from each GPE with all other objects
        NewVals[i+2] and NewVals[i+3] sends the x and y positions for the ith and jth object into the accel function to get a return of the new acceleration for the ith object
        NewVals[i] and NewVals[i+1] get the x and y velocities passed into them
        '''
        for i in range (0,NumofObjects * 4 -1,4):
    
            for j in range(0,NumofObjects * 4 - 1,4):
                if i != j:
                
                    NewVals[i+2] += Accel(vals[i],vals[j], vals[i+1], vals[j+1], Masses[j//4])
                    NewVals[i+3] += Accel(vals[i+1],vals[j+1], vals[i], vals[j], Masses[j//4])
                    
       
            NewVals[i] = vals[i+2]
            NewVals[i+1] = vals[i+3]
        
       
    
        return NewVals
    
    start_time = time.time()
    
    '''
    All the postiitons and velocities at each time step are stored in the 2D array FinalVals
    '''
    FinalVals = odeint(Func,Invals,t, rtol = 1e-13, atol = 1e-13)
    
    end_time = time.time()
    print('Time taken = ', end_time - start_time)
    
    axis = fig.add_subplot(4,2,ma)
    colours = np.array(['y*','k','m','g','b','r','c--'])
    
    
    '''
    loops in steps of 4 so that the ith value is the x pos and the i + 1 value is the y pos
    i//4 is which object is being plotted so it knows what colour to plot it as
    '''
    for i in range(0,NumofObjects * 4 -1, 4):
        
        axis.plot((FinalVals[:,i]- FinalVals[:,0]) / AU,(FinalVals[:,i+1] - FinalVals[:,1] )/ AU ,colours[i//4])
       # plt.axis('scaled')
        axis.set_title(Titles[ma-1])
        axis.set_xlabel('x [AU]')
        axis.set_ylabel('y [AU]')
    
    '''
    only plots the legend on the first plot
    fancybox and shadow are aestetics
    bbox_to_anchor moves the legend off of the plot and into the middle of the figure
    '''
    
    if ma == 1:
        plt.legend(('Sun','Jupiter','Saturn','Uranus','Neptune','Pluto','New Star'), fancybox = True, shadow = True, bbox_to_anchor=(1,1.25))

    '''
    calculating the change in energy, each object needs the KE but only need to work out the GPE for each pair once
    '''
    TotalE = 0
    
    for i in range (0,NumofObjects * 4 -1,4):
        V = np.sqrt((FinalVals[:,i + 2])**2 + (FinalVals[:,i + 3])**2)
        TotalE += 0.5 * Masses[i//4] * V**2
        for j in range(i+4,NumofObjects * 4 - 3,4):
                R  = np.sqrt((FinalVals[:,i] - FinalVals[:,j])**2 + (FinalVals[:,i+1] - FinalVals[:,j+1])**2) 
                TotalE += -G * Masses[i//4] * Masses[j//4] / R
                
    deltaE = (TotalE - TotalE[0] ) / TotalE[0]

    
    
    Fig3 = plt.figure()
    plt.plot(t/secsinyear, deltaE)

