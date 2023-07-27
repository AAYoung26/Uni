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

start_time = time.time()
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
setting up the initial mass of the sun as well as creating an array to store the averge ditances from it
'''



ii = 10
RSun = np.zeros(shape =(ii,6))
Mi = 0.2
colours = np.array(['r','m','g','b','y--','c--'])
for i in range(1,ii):
    total = 0


    '''
    creating a class that can make any of the large objects like a planet/sun
    '''
    
    class Star:
        def __init__(self,Mass,Luminosity,Multiplier):
            self.x = 0
            self.y = 0
            self.Vx = 0
            self.Vy = 0
            self.mass = Mass * M0 * Multiplier
            self.Lum = L0 * (Luminosity * Multiplier)**3.5
            
    Sun = Star(Mi,Mi,i)
    
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
    Jupiter = SpaceObj(5.2,4331,1.9e27,0.34)
    Saturn = SpaceObj(9.5,10747,5.68e26,0.34)
    
    MyObj = np.array([Sun,Mercury,Venus,Earth,Mars])#,Jupiter,Saturn]) #storing all of the objects in one array
    
    
    '''
    Creating a function to work out and return the boundaries for the goldilocks zone
    '''

    def GoldilocksZone(Star):
        iZone = np.sqrt((Star.Lum/L0)/1.1) 
        oZone = np.sqrt((Star.Lum/L0)/0.53) 
        return iZone,oZone


    '''
    making the boundaries in m instead of AU
    '''

    iZ,Oz = GoldilocksZone(Sun) 
    iZ *= AU
    Oz *= AU


    '''
    creating arrays that will store the x,y,Vx,Vy vals of all the objects, also creating one for the masses
    '''
    Invals = np.array([])
    Masses = np.array([])
    
    NumofObjects = len(MyObj)
    for j in range(0,len(MyObj)):
        Invals = np.append(Invals,MyObj[j].x)
        Invals = np.append(Invals,MyObj[j].y)
        Invals = np.append(Invals,MyObj[j].Vx)
        Invals = np.append(Invals,MyObj[j].Vy)
        Masses = np.append(Masses,MyObj[j].mass)
    
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
       
        for j in range (0,NumofObjects * 4 -1,4):
            B = np.arange(0,NumofObjects*4,4, dtype = int)
            B = np.delete(B,int(j/4))
            iarr = np.ones(shape = np.shape(B), dtype = int ) * j
    
            NewValsx = Accel(vals[iarr],vals[B], vals[iarr+1], vals[B+1], Masses[B//4]) 
            NewValsy = Accel(vals[iarr+1],vals[B+1], vals[iarr], vals[B], Masses[B//4])
            
            NewVals[j+2] = np.sum(NewValsx, axis = 0)
            NewVals[j+3] = np.sum(NewValsy, axis = 0)
       
            NewVals[j] = vals[j+2]
            NewVals[j+1] = vals[j+3]
               
        return NewVals
    

     
    '''
    All the postiitons and velocities at each time step are stored in the 2D array FinalVals
    '''
    FinalVals = odeint(Func,Invals,t, rtol = 1e-12, atol = 1e-12)
    
   
   
    
    
    
    

    
    '''
    loops in steps of 4 so that the j + 1 value is the x pos and the j + 2 value is the y pos, it is done this way so that the Sun's avg distance from itself doesnt get calculated'
    j//4 is which object is being plotted so it knows which planet it is finding the avg ditance for 
    '''
    for j in range(0,NumofObjects * 4 -1, 4):
        
       RSun[i-1,j//4] = np.sum(np.sqrt((FinalVals[:,j+1]- FinalVals[:,0])**2 + (FinalVals[:,j+2] - FinalVals[:,1] )**2),axis = 0) / 183
        
    RSun[i-1,4] = iZ
    RSun[i-1,5] = Oz
for j in range(0,6):
    plt.plot(np.arange(0,ii) * Mi,RSun[:,j], colours[j])
'''
fancybox and shadow are aestetics
'''    
plt.legend(('Mercury', 'Venus', ' Earth','Mars','inner radius of Goldilocks zone','Outer radius of Goldilocks zone'),fancybox = True, shadow = True)
plt.xlabel('Mass of Sun in solar Masses')
plt.ylabel('Average distance from the Sun (m)')
plt.title('The average distance from the Sun with varying masses and the Goldilocks Zone')


end_time = time.time()
print('Time taken = ', end_time - start_time)
