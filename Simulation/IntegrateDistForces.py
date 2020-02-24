# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 18:44:03 2020

@author: Burhan
"""


import numpy as np
import GlobalConstants as g
import aerodynamicloads as loads
from matplotlib import pyplot as plt

A, B, Force_Distrib, Torque_Distrib = loads.DistForceTorqueMatrix()

l = len(Force_Distrib)

CFs = np.zeros((5, l, 8))
CTs = np.zeros((5, l, 8))

CFs[0, :, :4]  = Force_Distrib[:, :4]
CTs[0, :, :4]  = Torque_Distrib[:, :4]

xnodes = Force_Distrib[:,4]

def setcoefs(Cs):

    '''This function sets the coefficients for different levels of integration of the 
       force/torque distributions'''
  
    for i in range (l):
        
        for j in range (4):
            
            e = 3 - j
            
            for k in range(1, 5):
            
                Cs[k, i, j] = Cs[k - 1, i, j]/(e+k)

        if (i>0):
            
            for k in range(4+1):
                
                Cs[1, i, 4] += ((xnodes[i] - xnodes[i-1])**(4-k))*(Cs[1, i-1, k])
                
            Cs[2, i, 4] = Cs[1, i, 4]
            Cs[3, i, 4] = Cs[1, i, 4]/2
            Cs[4, i, 4] = Cs[1, i, 4]/6
            
        if (i>1):
            
            for k in range(4+2):
                
                Cs[2, i, 5] += ((xnodes[i] - xnodes[i-1])**(5-k))*(Cs[2, i-1, k])
                
            Cs[3, i, 5] = Cs[2, i, 5]
            Cs[4, i, 5] = Cs[2, i,5]/2
            
        if (i>2):
            
            for k in range(4+3):
                
                Cs[3, i, 6] += ((xnodes[i] - xnodes[i-1])**(6-k))*(Cs[3, i-1, k])
                
            Cs[4, i, 6] = Cs[3, i, 6]
            
        if (i>3):
            
            for k in range(4+4):
                
                Cs[4, i, 7] += ((xnodes[i] - xnodes[i-1])**(7-k))*(Cs[4, i-1, k])
  

      
def setallcoefs():
    
    A, B, Force_Distrib, Torque_Distrib = loads.DistForceTorqueMatrix()

    l = len(Force_Distrib)

    global CFs = np.zeros((5, l, 8))
    global CTs = np.zeros((5, l, 8))

    global Fs[0, :, :4]  = Force_Distrib[:, :4]
    global CTs[0, :, :4]  = Torque_Distrib[:, :4]

    global xnodes = Force_Distrib[:,4]
    
    setcoefs(CFs)
    setcoefs(CTs)




def w(x, il):
    
    '''This function calculates the value of the (il)^th integral of the 
       force distribution at x'''
    
    w = 0
    
    for i in range (1,41):
    
        if xnodes[i-1] < x <= xnodes[i]:
            
            for j in range(4+il):
                
                e = 3 - j + il
                
                w += CFs[il, i-1, j]*((x - xnodes[i-1])**e)
                
            break
    
    return(w)

def t(x, il):
    
    '''This function calculates the value of the (il)^th integral of the 
       torque distribution at x'''
    
    w = 0
    
    for i in range (1,41):
    
        if xnodes[i-1] < x <= xnodes[i]:
            
            for j in range(4+il):
                
                e = 3 - j + il
                
                w += CTs[il, i-1, j]*((x - xnodes[i-1])**e)
                
            break
    
    return(w)









'''
x = np.linspace(xnodes[0], xnodes[-1], 200)
dx = x[1]-x[0]

y0 = np.zeros(len(x))
y1 = np.zeros(len(x))
y2 = np.zeros(len(x))
y3 = np.zeros(len(x))
y4 = np.zeros(len(x))

for i in range (len(x)):
    y0[i] = t(x[i], 0)
    y1[i] = t(x[i], 1)
    y2[i] = t(x[i], 2)
    y3[i] = t(x[i], 3)
    y4[i] = t(x[i], 4)
    
S1 = 0
S2 = 0
S3 = 0
S4 = 0
for j in range(len(x)-1):
    
    S1 += (y0[j] + y0[j+1])*0.5*dx
    S2 += (y1[j] + y1[j+1])*0.5*dx
    S3 += (y2[j] + y2[j+1])*0.5*dx
    S4 += (y3[j] + y3[j+1])*0.5*dx
    
print(S1, y1[-1])
print(S2, y2[-1])
print(S3, y3[-1])
print(S4, y4[-1])
'''
    
plt.plot(x,y0)
plt.plot(x,y1)
plt.plot(x,y2)
plt.plot(x,y3)
plt.plot(x,y4)
plt.grid(True)
plt.show()

            
            
    
        
        


