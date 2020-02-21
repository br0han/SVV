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
    
    if x >= Force_Distrib[41, 4]:
        
        x = Force_Distrib[41, 4]
        
        Coef = Force_Distrib[40, :4]
        
        for j in range(4):
            
            e = 3 - j
                
            w += (Coef[j]*(x**(e + int_level)))/(ncr(e+int_level, e))
            
    elif x <= Force_Distrib[0, 4]:
        
        x = Force_Distrib[0, 4]
        
        Coef = Force_Distrib[0, :4]
        
        for j in range(4):
                
            e = 3 - j
                
            w += (Coef[j]*(x**(e + int_level)))/(ncr(e+int_level, e))
                
    else:
        for i in range(1,41):
            
            if Force_Distrib[i-1, 4] < x < Force_Distrib[i, 4]:
                
                Coef = Force_Distrib[i-1, :4]
                
                for j in range(4):
                    
                    e = 3 - j
                
                    w += (Coef[j]*(x**(e + int_level)))/(ncr(e+int_level, e))
                
    return (w)

x = np.linspace(0, g.la, 1000)

y = W(x, 0)
    
plt.plot(x, y)
plt.show()
    
    
