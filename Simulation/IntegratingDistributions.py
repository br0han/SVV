# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:54:49 2020

@author: Burhan

This file contains the equilibirium equations?
"""

import numpy as np
import GlobalConstants as g
import aerodynamicloads as loads
from matplotlib import pyplot as plt

import operator as op
from functools import reduce

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer / denom


#import the force interpolant


"""Internal distributions"""

A, B, Force_Distrib, Torque_Distrib = loads.DistForceTorqueMatrix()


def W(x, int_level):
    '''Takes the matrix of the splines of the force distribution and converts
       them into one function of x
    '''
    
    w = 0
    
    '''if x >= Force_Distrib[40-1, 4]:
        
        x = Force_Distrib[40-1, 4]
        
        Coef = Force_Distrib[40-1, :4]
        
        for j in range(4):
            
            e = 3 - j
                
            w += (Coef[j]*(x**(e + int_level)))/(ncr(e+int_level, e))
            
            
    elif x <= Force_Distrib[0, 4]:
        
        x = Force_Distrib[0, 4]
        
        Coef = Force_Distrib[0, :4]
        
        for j in range(4):
                
            e = 3 - j
                
            w += (Coef[j]*(x**(e + int_level)))/(ncr(e+int_level, e))'''
                
    #else:
    
    C = 0
    
    for i in range(1,40):
            
        if Force_Distrib[i-1, 4] < x <= Force_Distrib[i, 4]:
            
            x0 = Force_Distrib[i-1, 4]
                
            Coef = Force_Distrib[i-1, :4]
            
            for j in range(4):
                
                e = 3 - j
            
                w += (Coef[j]*((x - x0)**(e + int_level)))/(ncr(e+int_level, e))
                
            break
                
        if (i>1 and int_level>0):
            
            x0prev = Force_Distrib[i-1, 4]
            x0 = Force_Distrib[i, 4]
            Coef_prev = Force_Distrib[i-1, :4]
            Ctemp = 0
            
            for j in range(4):
                
                e = 3 - j
                Ctemp += (Coef_prev[j]*((x0 - x0prev)**(e + int_level)))/(ncr(e+int_level, e))
                
            C+=Ctemp
                
    return (w + C)


    '''redoing the whole thing taking x to be an array'''
    
    #w = np.array(len(x))
    
    #for i in range (len(x)):
        
    
    
    
'''---------------------------------------------------------------------------------------------'''

x = np.linspace(Force_Distrib[1,4], Force_Distrib[39,4], 10000)
y0 = np.zeros(len(x))
y1 = np.zeros(len(x))
y2 = np.zeros(len(x))
y3 = np.zeros(len(x))
y4 = np.zeros(len(x))

for i in range(len(x)):
    y0[i] = W(x[i], 0)
    y1[i] = W(x[i], 1)
    y2[i] = 15*W(x[i], 2)
    y3[i] = 100*W(x[i], 3)
    y4[i] = 500*W(x[i], 4)
    
    
plt.plot(x, y0)
plt.plot(x, y1)
plt.plot(x, y2)
plt.plot(x, y3)
plt.plot(x, y4)

plt.legend((y1, y2, y3, y4), ('1','2','3','4'))

plt.grid(1)
plt.show()
    
    
