# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 14:57:01 2020

@author: Burhan
"""



import numpy as np
import GlobalConstants as g
import aerodynamicloads as loads
from matplotlib import pyplot as plt


class DistChordIntegrator:
    
    l = 83
    
    def __init__(self, C0cs, nodes):
        
        self.l = 83
        self.Ccs = np.zeros((5, self.l - 1, 10))
        self.Ccs[0, :, :4] = C0cs
        self.nodes = nodes #np.zeros((l, 1)) 

    def setcoefs(self, f = 0):
    
        '''This function sets the coefficients for different levels of integration of the 
           force/torque distributions'''
      
        for i in range (self.l - 1):
            
            for j in range (4):
                
                e = 3 - j + f
                
                for k in range(1, 5):
                
                    self.Ccs[k, i, j] = self.Ccs[k - 1, i, j]/(e+k)
    
            if (i>0):
                
                for k in range(4+1):
                    
                    self.Ccs[1, i, 4] += ((self.nodes[i] - self.nodes[i-1])**(4-k + f))*(self.Ccs[1, i-1, k])
                    
                self.Ccs[2, i, 4] = self.Ccs[1, i, 4]
                self.Ccs[3, i, 4] = self.Ccs[1, i, 4]/2
                self.Ccs[4, i, 4] = self.Ccs[1, i, 4]/6
                
            if (i>1):
                
                for k in range(4+2):
                    
                    self.Ccs[2, i, 5] += ((self.nodes[i] - self.nodes[i-1])**(5-k))*(self.Ccs[2, i-1, k])
                    
                self.Ccs[3, i, 5] = self.Ccs[2, i, 5]
                self.Ccs[4, i, 5] = self.Ccs[2, i,5]/2
                
            if (i>2):
                
                for k in range(4+3):
                    
                    self.Ccs[3, i, 6] += ((self.nodes[i] - self.nodes[i-1])**(6-k))*(self.Ccs[3, i-1, k])
                    
                self.Ccs[4, i, 6] = self.Ccs[3, i, 6]
                
            if (i>3):
                
                for k in range(4+4):
                    
                    self.Ccs[4, i, 7] += ((self.nodes[i] - self.nodes[i-1])**(7-k))*(self.Ccs[4, i-1, k])
                    
    def setallcoefs(self, f = 0):
        self.setcoefs()
    
    def w(self, x, il, f = 0):
        
        '''This function calculates the value of the (il)^th integral of the 
           force distribution at x'''
        
        w = 0
        
        for i in range (1,self.l):
        
            if self.nodes[i-1] < x <= self.nodes[i]:
                
                for j in range(4+il):
                    
                    e = 3 - j + il + f
                    
                    
                    w += self.Ccs[il, i-1, j]*((x - self.nodes[i-1])**e)
                   
                    
                break
        
        return(w)
    
#c1 = np.ones((82, 4))
    
    
#integ = DistChordIntegrator(c1)
    
    