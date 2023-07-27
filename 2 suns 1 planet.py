# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 16:51:55 2022

@author: ppyay1
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 13:02:45 2022

@author: Aaron
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time 



G = 6.67 * 10**(-11) #constant can be changed to 6.67x10^-11 for realism
secsinyear = 3.154*10**7

t = np.arange(0,20 * secsinyear, secsinyear/(365))

AU = 1.496e11

M0 = 1.989e30
ly = 9.46e15
Day1 = 60*60*24

star_x = -0.5 * AU
star_y  = 0 * AU
SMass1 = 1 * M0

Star2_x = 0.5 * AU
Star2_y = 0 * AU
SMass2 = 1 * M0

Planet_x = 0
Planet_y = 1 * AU
PMass = 5.97e24



Masses = np.array([SMass1,SMass2,PMass])


NumofObjects = len(Masses)

Invals = np.array([star_x,star_y,0, -15000,Star2_x,Star2_y,0,15000, Planet_x,Planet_y,40000,0])

def Accel(Pos1x,Pos2x,Pos1y,Pos2y,Mass):
    r1 = Pos1x - Pos2x
    r2 = Pos1y - Pos2y

    R = np.sqrt(r1**2 + r2**2)
    Acc = -G * Mass * r1 / R**3

    return Acc

def Func(vals,dt):

    NewVals = np.zeros(NumofObjects * 4)
   
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
        
        
    
    for i in range (0,NumofObjects * 4 -1,4):

        for j in range(0,NumofObjects * 4 - 1,4):
            if i != j:
                k = int(j/4)
                NewVals[i+2] += Accel(vals[i],vals[j], vals[i+1], vals[j+1], Masses[k])
                NewVals[i+3] += Accel(vals[i+1],vals[j+1], vals[i], vals[j], Masses[k])
                
   
        NewVals[i] = vals[i+2]
        NewVals[i+1] = vals[i+3]
    
   

    return NewVals
start_time = time.time()
FinalVals = odeint(Func,Invals,t, rtol = 1e-10, atol = 1e-10)

end_time = time.time()
print('Time taken = ', end_time - start_time)

Fig1 = plt.figure()
plt.subplots_adjust(wspace=1, hspace=0.5, top=0.85, left = 0.15)
plt.suptitle('Small body orbiting 2 Large bodies')
ax1 = Fig1.add_subplot(2,1,1)
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
ax2 = Fig1.add_subplot(2,1,2)

#for i in range(0,10949):
colours = np.array(['r','b','g'])
for i in range(0,NumofObjects * 4 -1, 4):
    
    ax1.plot((FinalVals[:,i]) / AU,(FinalVals[:,i+1] )/ AU ,colours[i//4])
    ax1.axis('scaled')
    ax1.legend(('Body 1','Body 2', 'Body 3'))
# plt.plot((FinalVals[:,4] - FinalVals[:,0]) / AU,(FinalVals[:,5]-FinalVals[:,1]) / AU ,'b')  # can subtract by the position of the first star so that it is in a fixed position
# plt.plot((FinalVals[:,8]-FinalVals[:,0]) / AU,(FinalVals[:,9]-FinalVals[:,1]) / AU ,'g')
# plt.plot((FinalVals[:,12]-FinalVals[:,0]) / AU,(FinalVals[:,13]-FinalVals[:,1]) / AU ,'p')
# plt.plot((FinalVals[:,16]-FinalVals[:,0]) / AU,(FinalVals[:,17]-FinalVals[:,1]) / AU ,'r*')
# plt.plot(FinalVals[:,0] / AU,FinalVals[:,1] / AU ,'y*')
# plt.plot(FinalVals[:,4] / AU,FinalVals[:,5] / AU ,'b')  
# plt.plot(FinalVals[:,8]/ AU,FinalVals[:,9] / AU ,'g')
# plt.plot(FinalVals[:,12] / AU,FinalVals[:,13] / AU ,'p')
# plt.plot(FinalVals[:,16] / AU,FinalVals[:,17] / AU ,'r*')
# #plt.pause(0.0001)
#plt.cla()






ax2.plot(t/secsinyear,FinalVals[:,0]/AU,'r')
ax2.plot(t/secsinyear,FinalVals[:,4]/AU,'b')
ax2.plot(t/secsinyear,FinalVals[:,8]/AU,'g')

ax2.axis('scaled')
plt.xlabel('t [yr]')
plt.ylabel('x [AU]')
ax2.legend(('Body 1','Body 2', 'Body 3'))
TotalE = 0

for i in range (0,NumofObjects * 4 -1,4):
    V = np.sqrt((FinalVals[:,i + 2])**2 + (FinalVals[:,i + 3])**2)
    m = int(i/4)
    TotalE += 0.5 * Masses[m] * V**2
    for j in range(i+4,NumofObjects * 4 - 3,4):
            R  = np.sqrt((FinalVals[:,i] - FinalVals[:,j])**2 + (FinalVals[:,i+1] - FinalVals[:,j+1])**2) 
            n=int(j/4)
            TotalE += -G * Masses[m] * Masses[n] / R
            
deltaE = (TotalE - TotalE[0] ) / TotalE[0]

Fig3 = plt.figure()
plt.plot(t/secsinyear, deltaE)


# Sv = np.sqrt((FinalVals[:,2])**2 + (FinalVals[:,3])**2)    
    
# S1KE = 0.5 * SMass1 * (Sv)**2

# Pv = np.sqrt((FinalVals[:,6])**2 + (FinalVals[:,7])**2)    

# P1KE = 0.5 * PMass1 * Pv**2

# P2v = np.sqrt((FinalVals[:,10])**2 + (FinalVals[:,11])**2)    

# P2KE = 0.5 * PMass2 * P2v**2

# Av = np.sqrt((FinalVals[:,14])**2 + (FinalVals[:,15])**2)    

# AKE = 0.5 * AMass * Av**2

# TotalE = S1KE + P1KE + P2KE + AKE

# R  = np.sqrt((FinalVals[:,0] - FinalVals[:,4])**2 + (FinalVals[:,1] - FinalVals[:,5])**2)    

# S1P1 = -G *SMass1 * PMass1 / R

# R2  = np.sqrt((FinalVals[:,0] - FinalVals[:,8])**2 + (FinalVals[:,1] - FinalVals[:,9])**2)    

# S1P2 = -G *SMass1 * PMass2 / R2

# R3  = np.sqrt((FinalVals[:,0] - FinalVals[:,12])**2 + (FinalVals[:,1] - FinalVals[:,13])**2)    

# S1A = -G *SMass1 * AMass / R3

# TotalE += S1P1 + S1P2 + S1A

# R4  = np.sqrt((FinalVals[:,4] - FinalVals[:,8])**2 + (FinalVals[:,5] - FinalVals[:,9])**2)    

# P1P2 = -G *PMass1 * PMass2 / R4

# R5  = np.sqrt((FinalVals[:,4] - FinalVals[:,12])**2 + (FinalVals[:,5] - FinalVals[:,13])**2)    

# P1A = -G *PMass1 * AMass / R5

# TotalE += P1P2 + P1A

# R6  = np.sqrt((FinalVals[:,8] - FinalVals[:,12])**2 + (FinalVals[:,9] - FinalVals[:,13])**2)    

# P2A = -G *PMass2 * AMass / R6

# TotalE += P2A


