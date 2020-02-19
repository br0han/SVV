# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:53:42 2020

@author: Burhan

Methods to calculate the moments of inertia of the aileron profile
"""

import GlobalConstants as g
import numpy as np
import Coordinates as coors
import math

def I_ring():
    
    Izz = 0.5*np.pi*(g.ha/2.)^3
    
    Iyy = g.tsk*((g.ha/2)**3)*((np.pi**2 - 8)/(2*np.pi)) + (np.pi*(0.5*g.ha)*g.tsk)*((2*(0.5*g.ha)/np.pi) + centroid)**2
    
    return [Iyy, Izz]


def I_spar():
    
    Izz = g.tsp*(g.ha^3)/12
    
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
    
    for i in range(0, g.nst/2):
        
        Izz += 2*Ast*(ycoords[i])**2
        
    #Calculate Iyy 
    
    for j in range(0, g.nst):
        
        Iyy += Ast*((zcoords[i] + centroid)**2)
        
    return [Iyy, Izz]
        
    
def I_skin():
    
    alpha = math.atan((g.ha/2)/(g.Ca - 0.5*g.ha)) #angle of staight skin to z-axis
    
    lsk = math.sqrt((g.ha/2.)**2 + (g.Ca - 0.5*g.ha)**2) #length of straigth skin 
    
    Izz = 2*((g.tsk*(lsk**3)*((math.sin(alpha))**2))/12 + ((g.ha/4.)**2)*(lsk*g.tsk))
    
    Iyy = 2*((g.tsk*(lsk**3)*((math.cos(alpha))**2))/12 + (Ast*(-(0.5*(g.Ca - 0.5*g.ha)) + centroid)**2))
    
    return [Iyy, Izz]
    
    
def I_zz():
    
    return (I_ring[1] + I_spar[1] + I_str[1] + I_skin[1])


def I_yy():
    
    return (I_ring[0] + I_spar[0] + I_str[0] + I_skin[0])
    
    
    
    