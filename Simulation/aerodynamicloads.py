import GlobalConstants as GC
import numpy as np
from cubicspline import cubicspline
import MyNewIntegralForTorque as icf
import Shear_center as sc

import time


def DistForceTorqueMatrix():
    
    print ("Called")
    
    data = np.genfromtxt('aerodynamicloadcrj700.dat', dtype=None, delimiter=',')
    
    ShearCenter = sc.FindSC()  ##changed to 0 from 0.007
    
    print("SC", ShearCenter)
    
    Ca = GC.Ca #0.484
    la = GC.la #1.691
    ha = GC.ha
    
    Coord = []
    for i in range(1,82):
            theta0 = (i-1)/81 * np.pi
            theta1 = (i)/81 * np.pi
            z = -1/2*(Ca/2*(1-np.cos(theta0))+Ca/2*(1-np.cos(theta1)))
            Coord.append(-z)
    Coord = np.array(Coord)
    Coord = np.append(0,Coord)
    Coord = np.append(Coord,Ca)
        ###End of inputs
    
    Span = []
    for i in range(1,42):
            theta0 = (i-1)/41 * np.pi
            theta1 = (i)/41 * np.pi
            z = 1/2*(la/2*(1-np.cos(theta0))+la/2*(1-np.cos(theta1)))
            Span.append(z)
    Span = np.array(Span)
    Span = np.append(0,Span)
    Span = np.append(Span,la)
    
    Force_Distrib = np.zeros(41)
    Torque_Distrib = np.zeros(41)
    
    for i in range (0, 41):
        val = data[:,i]
        val = np.append(0,val)
        val = np.append(val,0)
        
        CoefM = cubicspline(val, Coord)
        
        integ = icf.DistChordIntegrator(CoefM, Coord)
        integ.setallcoefs()
        
        #Force calc
        Force_Distrib[i] = integ.w(Ca, 1)
        
        #Torque calc
        x = np.linspace(0, Ca, 3000)
    
        dx = Ca/len(x)
   
        F0 = np.zeros(len(x))
        
        for k in range (len(x)):
            
            F0[k] = integ.w(x[k], 0)
        
        Fx = 0
        for j in range (len(x) - 1):
            Fx += ((F0[j] + F0[j+1])*0.5)*(x[j]+0.5*dx)*dx
       
                   
        PoA = -(Fx/integ.w(Ca,1) - ha/2)
        
        D2sc = PoA - ShearCenter
        
        Torque_Distrib[i] = integ.w(Ca, 1)*D2sc
           
    Force_Distrib = np.append(0, Force_Distrib)
    Force_Distrib = np.append(Force_Distrib, 0)
    
    Torque_Distrib = np.append(0, Torque_Distrib)
    Torque_Distrib = np.append(Torque_Distrib, 0)
    
    FD = cubicspline(Force_Distrib, Span)
    TD = cubicspline(Torque_Distrib, Span)
    
    Span  = np.reshape(Span,(len(Span),1))
    BotRow = np.array([0,0,0,0])
    
    FD = np.vstack((FD,BotRow))
    TD = np.vstack((TD,BotRow))
    
    FD = np.hstack((FD,Span))
    TD = np.hstack((TD,Span))  

    
    return (FD, TD)

