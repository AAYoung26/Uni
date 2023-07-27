# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 16:24:58 2022

@author: Aaron
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

fig = plt.figure()
#plt.subplots_adjust(wspace=0.6, hspace=1, top=0.85, left = 0.15)
plt.suptitle('Moon orbit for 1 year when a large object is nearby')
neededJets = 0
Found = False
JetVals = np.array([0,1e-05,3e-05,5e-05,7e-05,9e-05])
MassVals = np.array([0.0001,0.001,0.01,0.1,1])

G = 6.67 * 10**(-11) #constant can be changed to 6.67x10^-11 for realism
secsinyear = 3.154*10**7

t = np.arange(0,1 * secsinyear, secsinyear/(365*24*60))

AU = 1.496 * 10**11

M0 = 1.989*10**30
ly = 9.46 * 10**15
Day1 = 60*60*24

Sun_x = 0 * AU
Sun_y  = 0 * AU
SMass1 = 1 * M0

Earth_x = 0.98 * AU
Earth_y = 0 * AU
EMass = 5.97*10**24
EPeriod = Day1 * 365.2
EVel = Earth_x * 2*np.pi/EPeriod

Moon_x = Earth_x + 384400000
Moon_y = 0 * AU
MMass = 7.346*10**22
MPeriod = 27.3*Day1
MVel = 384400000 * 2 * np.pi / (MPeriod)

Ast_x = (Earth_x + 384400000*2)
Ast_y = -50000000
AMass = EMass * 0.001


#for i in range(0,len(MassVals)):
Masses = np.array([EMass,MMass, AMass])


NumofObjects = len(Masses)
    

 #   Jets = JetVals[i]
    
Invals = np.array([Earth_x,Earth_y,0,0,Moon_x,Moon_y,0,-100000, Ast_x,Ast_y,0,0])

def Accel(Pos1x,Pos2x,Pos1y,Pos2y,Mass):
    r1 = Pos1x - Pos2x
    r2 = Pos1y - Pos2y

    R = np.sqrt(r1**2 + r2**2)
    Acc = -G * Mass * r1 / R**3

    return Acc
    
      #  Jets += 0.0000001
        
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
        
    #drag = 0.0001175 
  
    
    # angle = tan-1 (dx/dy)
    
    for i in range (0,NumofObjects * 4 -1,4):

        for j in range(0,NumofObjects * 4 - 1,4):
            if i != j:
                k = int(j/4)
                NewVals[i+2] += Accel(vals[i],vals[j], vals[i+1], vals[j+1], Masses[k])
                NewVals[i+3] += Accel(vals[i+1],vals[j+1], vals[i], vals[j], Masses[k])
  
        NewVals[i] = vals[i+2]
        NewVals[i+1] = vals[i+3]
       # Acc = np.sqrt(NewVals[i+2] **2 + NewVals[i+3] **2)
  #  dx  =vals[0] - vals[4]
   # dy = vals[1] - vals[5]
  #  ang = np.arctan(dx/dy)
  #  Jetsx = Jets * np.sin(ang) 
  #  Jetsy = Jets * np.cos(ang)
    
  #  NewVals[6] += Jetsx
  #  NewVals[7] += Jetsy

        # if NewVals[i+2] >= 0:
        #     NewVals[i+2]-= drag
        # else:
        #     NewVals[i+2] += drag
            
        # if NewVals[i+3] >= 0:
        #     NewVals[i+3]-= drag
        # else:
        #     NewVals[i+3] += drag
   
            
    
        return NewVals
start_time = time.time()
FinalVals = odeint(Func,Invals,t, rtol = 1e-12, atol = 1e-12)
#R = np.sqrt((FinalVals[:,0]-FinalVals[:,4])**2 + (FinalVals[:,1]-FinalVals[:,5])**2)
# if np.min(R) < 0.001 * AU and np.max(R) < 0.1 * AU:
#         neededJet = Jets
#         Found = True
#         print(neededJet)
# if Jets > 1:
#     Found = True
#     print('no match')for i in range(0,NumofObjects * 4 -1, 4):
colours = np.array(['bo','k','r'])
#axis = fig.add_subplot(4,2,i+1)
  #  Force = Jets * MMass
    
   # print(round(Force,-15))
   # axis.set_title(round(Force,-15))
for i in range(0,NumofObjects * 4 -1, 4):
    plt.plot((FinalVals[:,i]-FinalVals[:,0]) / AU,(FinalVals[:,i+1]-FinalVals[:,1]) / AU ,colours[i//4])
end_time = time.time()
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
print('Time taken = ', end_time - start_time)
    

#Fig1 = plt.figure()

#for j in range(0,10949):
#colours = np.array(['bo','k'])
# for i in range(0,NumofObjects * 4 -1, 4):

#     plt.plot((FinalVals[:,i]-FinalVals[:,0]) / AU,(FinalVals[:,i+1]-FinalVals[:,1]) / AU ,colours[i//4])
    

#plt.axis('scaled')
#plt.xlabel('x [AU]')
#plt.ylabel('y [AU]')
plt.legend(('Earth', 'Moon', 'Large Object'))

TotalE = 0
Vel = np.sqrt((FinalVals[:,6])**2 + (FinalVals[:,7])**2)
SmallestVel = np.min(Vel)


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
