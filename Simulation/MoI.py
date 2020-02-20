# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:53:42 2020

@author: Burhan

Methods to calculate the moments of inertia of the aileron profile
"""

import GlobalConstants as g
import Centroid as cent
import numpy as np
import Coordinates as coors
import math

centroid = cent.Centroid_z()

ra = g.ha/2.

def I_ring():
    
    Izz = 0.5*np.pi*(ra)**3*g.tsk
    
    Iyy = g.tsk*((ra)**3)*((np.pi**2 - 8)/(2*np.pi)) + (np.pi*(ra)*g.tsk)*((2*(ra)/np.pi) - centroid)**2    
    
    return [Iyy, Izz]


def I_spar():
    
    Izz = g.tsp*(g.ha**3)/12
    
    Iyy = (g.ha*g.tsp)*(centroid**2)
    
    return [Iyy, Izz]
  
    
def I_str():
    
    Ast = (g.hst + g.wst)*g.tst
    
    Izz = 0
    Iyy = 0
    
    stcoords = coors.Coord_out()
    
    ycoords = stcoords[:, 1]
    zcoords = stcoords[:, 0]
    
    #Calculate Izz
    
    for i in range(0, int(g.nst/2)):
        
        Izz += 2*Ast*(ycoords[i])**2
        
    #Calculate Iyy 
    
    for j in range(0, g.nst):
        
        Iyy += Ast*((zcoords[j] - centroid)**2)
        
    return [Iyy, Izz]
        
    
def I_skin():
    
    alpha = math.atan(ra/(g.Ca - ra)) #angle of staight skin to z-axis
    
    lsk = math.sqrt(ra**2 + (g.Ca - ra)**2) #length of straigth skin 
    
    Izz = 2*((g.tsk*(lsk**3)*((math.sin(alpha))**2))/12 + ((ra/2.)**2)*(lsk*g.tsk))
    
    Iyy = 2*((g.tsk*(lsk**3)*((math.cos(alpha))**2))/12 + (lsk*g.tsk*(-(0.5*(g.Ca - ra)) - centroid)**2))
    
    return [Iyy, Izz]
    
    
def I_zz():
    
    return (I_ring()[1] + I_spar()[1] + I_str()[1] + I_skin()[1])


def I_yy():
    
    return (I_ring()[0] + I_spar()[0] + I_str()[0] + I_skin()[0])



print (I_yy(), I_zz())
    
    
    
    