# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:25:44 2020

@author: Burhan
"""
from math import sin,cos, radians
import numpy as np
import GlobalConstants as g
import MoI as moi
import Jcalc as JJ
import IntegrateDistForces as idf
import BigBoiMatrix as solvedunkowns
import Shear_center as SC

lol, Cs, FY, FZ = solvedunkowns.getunknowns()


Ry1 = Cs[0]; Rz1 = Cs[1]; Ry2 = Cs[2]; Rz2 = Cs[3]; Ry3 = Cs[4]; Rz3 = Cs[5]; Ra  = Cs[6]
Cv1 = Cs[7]; Cv2 = Cs[8]; Cw1 = Cs[9]; Cw2 = Cs[10]; Ct = Cs[11]

E = g.E
G = g.G

la = g.la; x1 = g.x1; x2 = g.x2; x3 = g.x3; xa = g.xa; ha = g.ha

d1 = g.d1
d3 = g.d3
P = g.Load2
theta_a = radians(g.maxdef)

# print ("Ry", Cs[0] + Cs[2] + Cs[4] + Cs[6]*sin(theta_a))
# print ("FY", FY)
# print ("Rz", Cs[1] + Cs[3] + Cs[5] + Cs[6]*cos(theta_a))
# print ("FZ", FZ)

sc = SC.FindSC() #get Shear centre from marianos code

Py = P*sin(theta_a)
Pz = P*cos(theta_a)

Ray = Ra*sin(theta_a)
Raz = Ra*cos(theta_a)

GJ = (1/(G*(JJ.J())))       #Flexural rigidity stuff
FRzz = (1/(E*moi.I_zz()))
FRyy = (1/(E*moi.I_yy()))

integ = idf.DistIntegrator()  


def step(x, p):
    ''' Definition of Macaulay Step function'''
    if x < 0:
        return 0
    elif x >= 0:
        return (x**p)
    
    
def Mz(x):
    '''Internal moment about z'''
    return (-Ry1*step(x - x1, 1) - Ray*step(x - (x2 - xa/2), 1) - Ry2*step(x - x2, 1) + Py*step(x - (x2 + xa/2), 1) - Ry3*step(x - x3, 1) + integ.w(x,2))

def My(x):
    '''Internal moment about y'''
    return (-Rz1*step(x - x1,1) - Raz*step(x - (x2 - xa/2), 1) - Rz2*step(x - x2, 1) + Pz*step(x - (x2 + xa/2), 1) - Rz3*step(x - x3, 1))

def T(x):
    '''Internal Torque'''
    return -(Ry1*sc*step(x - x1, 0) + Raz*(ha/2)*step(x - (x2 - xa/2), 0) - Ray*(ha/2 - sc)*step(x - (x2 - xa/2), 0) + Ry2*sc*step(x - x2, 0) - Pz*(ha/2)*step(x - (x2 + xa/2), 0) + Py*(ha/2 - sc)*step(x - (x2 + xa/2), 0) + Ry3*sc*step(x - x3, 0) + integ.t(x,1))

def Sy(x):
    '''Internal shear force, in y'''
    return (- Ry1*step(x - x1, 0) - Ray*step(x - (x2 - xa/2), 0) - Ry2*step(x - x2, 0) + Py*step(x - (x2 + xa/2), 0) - Ry3*step(x - x3, 0) + integ.w(x, 1))

def Sz(x):
    '''Internal shear force, in z'''
    return (-Rz1*step(x - x1,0) - Raz*step(x - (x2 - xa/2), 0) - Rz2*step(x - x2, 0) + Pz*step(x - (x2 + xa/2), 0) - Rz3*step(x - x3, 0))

def v(x):
    '''Deflection in y'''
    return 1000*(FRzz*((Ry1/6)*step(x - x1, 3) + (Ray/6)*step(x - (x2 - xa/2), 3) + (Ry2/6)*step(x - x2, 3) - (Py/6)*step(x - (x2 + xa/2), 3) + (Ry3/6)*step(x - x3, 3) - integ.w(x, 4) - Cv1*x - Cv2))

def w(x):
    '''Deflection in z'''
    return 1000*(FRyy*((Rz1/6)*step(x - x1, 3) + (Raz/6)*step(x - (x2 - xa/2), 3) + (Rz2/6)*step(x - x2, 3) - (Pz/6)*step(x - (x2 + xa/2), 3) + (Rz3/6)*step(x - x3, 3) - Cw1*x - Cw2))

def theta(x):
    '''Angular deflection'''
    return -(180/np.pi)*(GJ*(Ry1*sc*step(x - x1, 1) + Raz*(ha/2)*step(x - (x2 - xa/2), 1) - Ray*(ha/2 - sc)*step(x - (x2 - xa/2), 1) + Ry2*sc*step(x - x2, 1) - Pz*(ha/2)*step(x - (x2 + xa/2), 1) + Py*(ha/2 - sc)*step(x - (x2 + xa/2), 1) + Ry3*sc*step(x - x3, 1) + integ.t(x, 2) + Ct))


# print(v(x1)/cos(theta_a)*100, 100*v(x2), v(x3)/cos(theta_a)*100)

# print(w(x1)/sin(theta_a)*100, 100*w(x2), w(x3)/sin(theta_a)*100)


    
    
    
    