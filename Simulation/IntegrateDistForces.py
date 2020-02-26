# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 21:50:00 2020

@author: Burhan
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 18:44:03 2020

@author: Burhan
"""


import numpy as np
import GlobalConstants as g
import aerodynamicloads as loads
from matplotlib import pyplot as plt


class DistIntegrator:
    
    CFs = []
    CTs = []
    
    Force_Distrib, Torque_Distrib = loads.DistForceTorqueMatrix()
    
    l = len(Force_Distrib)
    
    CFs = np.zeros((5, l, 8))
    CTs = np.zeros((5, l, 8))
    
    CFs[0, :, :4]  = Force_Distrib[:, :4]
    CTs[0, :, :4]  = Torque_Distrib[:, :4]
    
    xnodes = Force_Distrib[:,4]

    def setcoefs(self, Cs):
    
        '''This function sets the coefficients for different levels of integration of the 
           force/torque distributions'''
      
        for i in range (self.l):
            
            for j in range (4):
                
                e = 3 - j
                
                for k in range(1, 5):
                
                    Cs[k, i, j] = Cs[k - 1, i, j]/(e+k)
    
            if (i>0):
                
                for k in range(4+1):
                    
                    Cs[1, i, 4] += ((self.xnodes[i] - self.xnodes[i-1])**(4-k))*(Cs[1, i-1, k])
                    
                Cs[2, i, 4] = Cs[1, i, 4]
                Cs[3, i, 4] = Cs[1, i, 4]/2
                Cs[4, i, 4] = Cs[1, i, 4]/6
                
            if (i>1):
                
                for k in range(4+2):
                    
                    Cs[2, i, 5] += ((self.xnodes[i] - self.xnodes[i-1])**(5-k))*(Cs[2, i-1, k])
                    
                Cs[3, i, 5] = Cs[2, i, 5]
                Cs[4, i, 5] = Cs[2, i,5]/2
                
            if (i>2):
                
                for k in range(4+3):
                    
                    Cs[3, i, 6] += ((self.xnodes[i] - self.xnodes[i-1])**(6-k))*(Cs[3, i-1, k])
                    
                Cs[4, i, 6] = Cs[3, i, 6]
                
            if (i>3):
                
                for k in range(4+4):
                    
                    Cs[4, i, 7] += ((self.xnodes[i] - self.xnodes[i-1])**(7-k))*(Cs[4, i-1, k])
                    
    def setallcoefs(self):
        self.setcoefs(self.CFs)
        self.setcoefs(self.CTs)

    def w(self, x, il):
        
        '''This function calculates the value of the (il)^th integral of the 
           force distribution at x'''
        
        w = 0
        
        for i in range (1,self.l):
        
            if self.xnodes[i-1] < x <= self.xnodes[i]:
                
                for j in range(4+il):
                    
                    e = 3 - j + il
                    
                    w += self.CFs[il, i-1, j]*((x - self.xnodes[i-1])**e)
                    
                break
        
        return(-w)
    
    def t(self, x, il):
        
        '''This function calculates the value of the (il)^th integral of the 
           torque distribution at x'''
        
        t = 0
        
        for i in range (1,self.l):
        
            if self.xnodes[i-1] < x <= self.xnodes[i]:
                
                for j in range(4+il):
                    
                    e = 3 - j + il
                    
                    t += self.CTs[il, i-1, j]*((x - self.xnodes[i-1])**e)
                    
                break
        
        return(t)



# integ = DistIntegrator()

# # integ.setcoefs(integ.CFs)

# integ.setallcoefs()


# x = np.linspace(integ.xnodes[0], integ.xnodes[-1], 10000)
# dx = x[1]-x[0]

# y0 = np.zeros(len(x))
# y1 = np.zeros(len(x))
# y2 = np.zeros(len(x))
# y3 = np.zeros(len(x))
# y4 = np.zeros(len(x))

# for i in range (len(x)):
#     y0[i] = integ.t(x[i], 0)
#     y1[i] = integ.t(x[i], 1)
#     y2[i] = integ.t(x[i], 2)
#     y3[i] = integ.t(x[i], 3)
#     y4[i] = integ.t(x[i], 4)
    
# S1 = 0
# S2 = 0
# S3 = 0
# S4 = 0
# for j in range(len(x)-1):
    
#     S1 += (y0[j] + y0[j+1])*0.5*dx
#     S2 += (y1[j] + y1[j+1])*0.5*dx
#     S3 += (y2[j] + y2[j+1])*0.5*dx
#     S4 += (y3[j] + y3[j+1])*0.5*dx
    
# print(S1, y1[-1])
# print(S2, y2[-1])
# print(S3, y3[-1])
# print(S4, y4[-1])

    
# plt.plot(x,y0)
# # plt.plot(x,y1)
# # plt.plot(x,y2)
# # plt.plot(x,y3)
# # plt.plot(x,y4)
# plt.grid(True)
# plt.show()

    
        
        


