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
plt.suptitle('Two heavy bodies in close orbit and a light body orbiting them both')

'''
defining some constants
'''
G = 6.67 * 10**(-11) 
secsinyear = 3.154*10**7

t = np.arange(0,15 * secsinyear, secsinyear/(365/4)) #creating an array of time values from 0 to 15 years in steps of a quarter of a year

AU = 1.496 * 10**11

M0 = 1.989*10**30

Day1 = 60*60*24


'''
creating a class that can make any of the large objects like a planet/sun
'''
class SpaceObj:
    def __init__(self,x,y,Vx,Vy,Mass):
        self.x = x * AU
        self.y = y * AU
        self.z = 0
        self.Vx = Vx
        self.Vy = Vy
        self.Vz = 0
        self.mass = Mass
    

Star1 = SpaceObj(-0.5,0,0,-15000,1 * M0)
Star2 = SpaceObj(0.5,0,0,15000,1*M0)
Planet = SpaceObj(0,1,40000,0,5.97e24)



MyObj = np.array([Star1,Star2,Planet]) #storing all of the objects in one array

'''
creating arrays that will store the x,y,z,Vx,Vy,Vz vals of all the objects, also creating one for the masses
'''
Invals = np.array([])
Masses = np.array([])

NumofObjects = len(MyObj)
for i in range(0,len(MyObj)):
    Invals = np.append(Invals,MyObj[i].x)
    Invals = np.append(Invals,MyObj[i].y)
    Invals = np.append(Invals,MyObj[i].z)
    Invals = np.append(Invals,MyObj[i].Vx)
    Invals = np.append(Invals,MyObj[i].Vy)
    Invals = np.append(Invals,MyObj[i].Vz)
    Masses = np.append(Masses,MyObj[i].mass)



'''
defining a function to work out the acceleration
'''
def Accel(Pos1x,Pos2x,Pos1y,Pos2y,Pos1z,Pos2z,Mass):
    r1 = Pos1x - Pos2x
    r2 = Pos1y - Pos2y
    r3 = Pos1z - Pos2z

    R = np.sqrt(r1**2 + r2**2 + r3**2)
    Acc = -G * Mass * r1 / R**3

    return Acc

'''
defining the function that will be used for the odeint cycle
'''

def Func(vals,dt):
    '''
    defining an empty array that will be filled with the new velocities and accelerations
    '''
    NewVals = np.zeros(NumofObjects * 6)
    
    '''
    B is an array of the indexs of the x pos for each object except for the ith object so that it doesnt do caluclations with itslef
    iarr is an array of the same shape as B but filled with the i value so that it is the index for the xpos of the ith object
    NewValsx and NewValsy send the x y and z positions for the ith and every other object into the accel function to get a return of an array of the acceleration contributions for the ith object
    NewVals[i+3] and NewVals[i+4] and NewVals[i+5] add the columns of NewValsx and NewValsy and NewValsz to get the total accel for the ith object
    NewVals[i] and NewVals[i+1] and NewVals[i+2] get the x and y and z velocities passed into them
    '''
   
    # for i in range (0,NumofObjects * 6 -1,6):
    #     B = np.arange(0,NumofObjects*6,6, dtype = int)
    #     B = np.delete(B,int(i/6))
    #     iarr = np.ones(shape = np.shape(B), dtype = int ) * i

    #     NewValsx = Accel(vals[iarr],vals[B], vals[iarr+1], vals[B+1],vals[iarr+2],vals[B+2], Masses[B//6])   #vectorised but it is slower for this number of bodies
    #     NewValsy = Accel(vals[iarr+1],vals[B+1], vals[iarr], vals[B],vals[iarr+2],vals[B+2], Masses[B//6])
    #     NewValsz = Accel(vals[iarr+2],vals[B+2], vals[iarr], vals[B],vals[iarr+1],vals[B+1], Masses[B//6])
        
    #     NewVals[i+3] = np.sum(NewValsx, axis = 0)
    #     NewVals[i+4] = np.sum(NewValsy, axis = 0)
    #     NewVals[i+5] = np.sum(NewValsz, axis = 0)
   
    #     NewVals[i] = vals[i+3]
    #     NewVals[i+1] = vals[i+4]
    #     NewVals[i+2] = vals[i+5]
        
    '''
    the outer loop gets the first object, the inner loop then works out the acceleration contribution from each GPE with all other objects
    NewVals[i+3] and NewVals[i+4] and NewVals[i+5] sends the x and y positions for the ith and jth object into the accel function to get a return of the new acceleration for the ith object
    NewVals[i] and NewVals[i+1] and NewVals[i+2] get the x and y and z velocities passed into them
    '''      
    
    for i in range (0,NumofObjects * 6 -1,6):

        for j in range(0,NumofObjects * 6 - 1,6):
            if i != j:
                k = int(j/6)
                NewVals[i+3] += Accel(vals[i],vals[j], vals[i+1], vals[j+1], vals[i+2],vals[j+2], Masses[k])
                NewVals[i+4] += Accel(vals[i+1],vals[j+1], vals[i], vals[j],vals[i+2],vals[j+2], Masses[k])
                NewVals[i+5] += Accel(vals[i+2],vals[j+2], vals[i], vals[j],vals[i+1],vals[j+1], Masses[k])
   
        NewVals[i] = vals[i+3]
        NewVals[i+1] = vals[i+4]
        NewVals[i+2] = vals[i+5]
   

    return NewVals
start_time = time.time()
'''
All the postiitons and velocities at each time step are stored in the 2D array FinalVals
'''
FinalVals = odeint(Func,Invals,t, rtol = 1e-12, atol = 1e-12)

end_time = time.time()
print('Time taken = ', end_time - start_time)

'''
setting up the axis for 3D plots
'''
ax = fig.add_subplot(121,projection ='3d')
colours = np.array(['b','r','g'])

'''
loops in steps of 6 so that the ith value is the x pos and the i + 1 value is the y pos
i//6 is which object is being plotted so it knows what colour to plot it as
'''
for i in range(0,NumofObjects * 6 -1, 6):
    
    ax.plot3D(FinalVals[:,i] / AU,FinalVals[:,i+1] / AU ,FinalVals[:,i+2] / AU , colours[i//6])

plt.axis('auto')
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
ax.set_zlabel('z [AU]')
plt.legend(('Star 1','Star 2', 'Planet'),fancybox = True, shadow = True)


'''
calculating the change in energy, each object needs the KE but only need to work out the GPE for each pair once
'''
TotalE = 0

for i in range (0,NumofObjects * 6 -1,6):
    V = np.sqrt((FinalVals[:,i + 3])**2 + (FinalVals[:,i + 4])**2 + (FinalVals[:,i+5]))
    m = int(i/6)
    TotalE += 0.5 * Masses[m] * V**2
    for j in range(i+6,NumofObjects * 6 - 5,6):
            R  = np.sqrt((FinalVals[:,i] - FinalVals[:,j])**2 + (FinalVals[:,i+1] - FinalVals[:,j+1])**2 + (FinalVals[:,i+2] - FinalVals[:,j+2])**2) 
            n=int(j/6)
            TotalE += -G * Masses[m] * Masses[n] / R
            
deltaE = (TotalE - TotalE[0] ) / TotalE[0]
ax2 = fig.add_subplot(122)
ax2.plot(t/secsinyear, deltaE)
ax2.set_xlabel('Time [yr]')
ax2.set_ylabel('$\Delta E$')

'''
calculatin the position from the x,y and z values and then plotting it against time for all 3 bodies
'''

Fig2 = plt.figure()
R = np.sqrt(FinalVals[:,0]**2 + FinalVals[:,1]**2 + FinalVals[:,2]**2)
plt.plot(t/secsinyear,R/AU,'b')

R = np.sqrt(FinalVals[:,6]**2 + FinalVals[:,7]**2 + FinalVals[:,8]**2)
plt.plot(t/secsinyear,R/AU,'r')

R = np.sqrt(FinalVals[:,12]**2 + FinalVals[:,13]**2 + FinalVals[:,14]**2)
plt.plot(t/secsinyear,R/AU,'g')

plt.axis('scaled')
plt.xlabel('Time [yr]')
plt.ylabel('R [AU]')
plt.legend(('Star 1','Star 2','Planet'),fancybox = True, shadow = True)
plt.title('Position of each object versus time')

