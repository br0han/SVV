# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:05:00 2020

@author: Burhan
"""

import numpy as np
import GlobalConstants as g
import IntegrateDistForces as idf
import MoI as moi
import Jcalc as JJ
from math import sin,cos,radians



def getunknowns():
    E = g.E
    G = g.G
    
    la = g.la
    x1 = g.x1
    x2 = g.x2
    x3 = g.x3
    xa = g.xa
    
    d1 = g.d1
    d3 = g.d3
    
    ha = g.ha
    
    P = g.Load2
    theta_a = radians(g.maxdef)
    
    sc = -0.0053559495 #get Shear centre from marianos code
    
    P_y = P*sin(theta_a)
    P_z = P*cos(theta_a)
    
    
    integ = idf.DistIntegrator()
    integ.setallcoefs()
    
    
    
    A = np.zeros((12,12))
    b = np.zeros((12, 1))
    
    
    
    '''Row 0 (Mz(la))'''
    A[0,0] = -(la - x1) #Ry1
    A[0,2] = -(la - x2) #Ry2
    A[0,4] = -(la - x3) #Ry3
    A[0,6] = -(la - (x2 - xa/2))*sin(theta_a) #Ray
    ''' RHS(0) '''
    b[0,0] = -P_y*(la - (x2 + xa/2)) - integ.w(la, 2)
    
    
    
    ''' Row 1 My(la)'''
    A[1,1] = -(la - x1) #Rz1
    A[1,3] = -(la - x2) #Rz2
    A[1,5] = -(la - x3) #Rz3
    A[1,6] = -(la - (x2 - xa/2))*cos(theta_a) #Raz
    '''RHS(1)'''
    b[1,0] = -P_z*(la - (x2 + xa/2))
    
    
    
    # '''Row 2 (T(la))'''
    # A[2,0] = sc #Ry1
    # A[2,2] = sc #Ry2
    # A[2,4] = sc #Ry3
    # A[2,6] = cos(theta_a)*(ha/2) - sin(theta_a)*(-sc + 0.5*ha) #Ra
    # '''RHS(2)'''
    # b[2,0] = -P_z*0.5*ha + P_y*(-sc + 0.5*ha) - integ.t(la,1)
    
    '''Row 2 (T(la)) = 0'''
    A[2,0] = sc #Ry1
    A[2,2] = sc #Ry2
    A[2,4] = sc #Ry3
    A[2,6] = cos(theta_a)*(ha/2) - sin(theta_a)*(-sc + 0.5*ha) #Ra
    '''RHS(2)'''
    b[2,0] = P_z*(ha/2) - P_y*(ha/2 - sc) - integ.t(la,1)
    
#made changes to Sy and Sz
    
    '''Row 3 (Sy(la)) = 0'''
    A[3,0] = -1 #Ry1
    A[3,2] = -1 #Ry2
    A[3,4] = -1 #Ry3
    A[3,6] = -sin(theta_a) #Ray
    '''RHS(3)'''
    b[3,0] = -P_y - integ.w(la, 1)
    
    
    
    '''Row 4 (Sz(la)) = 0'''
    A[4,1] = -1 #Rz1
    A[4,3] = -1 #Rz2
    A[4,5] = -1 #Rz3
    A[4,6] = -cos(theta_a) #Raz
    '''RHS(4)'''
    b[4,0] = -P_z 
    
    
    
    # '''Row 5 w(x1)'''
    # A[5,9] = la #Cw1
    # A[5,10] = 1 #Cw2
    # '''RHS(5)'''
    # b[5,0] = d1*sin(theta_a)
    
    '''Row 5 w(x1) = -d1sin(theta_a)'''
    FRyy = (1/(E*moi.I_yy()))
    
    ##A[5,1] = (FRyy/6)
    A[5,9] = -(FRyy)*x1 #Cw1
    A[5,10] = -(FRyy) #Cw2
    '''RHS(5)'''
    b[5,0] = -d1*sin(theta_a)
    
    
    
    # '''Row 6 w(x2)'''
    # FRyy = (1/(E*moi.I_yy()))
    # A[6,1] = -(1/6)*((x2 - x1)**3)*(-FRyy) #Rz1
    # A[6,6] = -(1/6)*((x2 - (x2 - 0.5*xa))**3)*(-FRyy)*cos(theta_a) #Raz
    # A[6,9] = x2 #Cw1
    # A[6,10] = 1 #Cw2
    # '''RHS(6)'''
    # #Zero
    
    '''Row 6 w(x2) = 0'''
    
    A[6,1] = (FRyy/6)*((x2 - x1)**3) #Rz1
    ##A[6,3] = (FRyy/6) #Ry2
    A[6,6] = -(FRyy/6)*((x2 - (x2 - 0.5*xa))**3)*cos(theta_a) #Raz
    A[6,9] = -(FRyy)*x2 #Cw1
    A[6,10] = -(FRyy) #Cw2
    '''RHS(6)'''
    #Zero
    
    
    
    # '''Row 7 w(x3) = d3sin(theta_a)'''
    # A[7,1] = -(1/6)*((x3 - x1)**3)*(-FRyy) #Rz1
    # A[7,3] = -(1/6)*((x3-x2)**3)*(-FRyy) #Rz2
    # A[7,6] = -(1/6)*((x3-(x2 - 0.5*xa))**3)*(-FRyy)*cos(theta_a) #Raz
    # A[7,9] = x3 #Cw1
    # A[7,10] = 1 #Cw2
    # '''RHS(7)'''
    # b[7,0] = -(P_z/6)*((x3 - (x2 + 0.5*xa))**3)*FRyy + d3*sin(theta_a)
    
    '''Row 7 w(x3) = -d3sin(theta_a)'''
    A[7,1] = (FRyy/6)*((x3 - x1)**3) #Rz1
    A[7,3] = (FRyy/6)*((x3 - x2)**3) #Rz2
    ##A[7,5] = (FRyy/6) #Rz3
    A[7,6] = (FRyy/6)*((x3 - (x2 - 0.5*xa))**3)*cos(theta_a) #Raz
    A[7,9] = -(FRyy)*x3 #Cw1
    A[7,10] = -(FRyy) #Cw2
    '''RHS(7)'''
    b[7,0] = (FRyy/6)*P_z*((x3 - (x2 + 0.5*xa))**3) - d3*sin(theta_a)
    
    
    
    # '''Row 8 v(x1) - sc*theta(x1)'''
    # GJ = (1/(G*(JJ.J())))
    # FRzz = (1/(E*moi.I_zz()))
    
    # A[8,7] = x1  #Cv1
    # A[8,8] = 1 #Cv2
    # A[8,11] = -sc #Ctheta
    # '''RHS(8)'''
    # b[8,0] = FRzz*integ.w(x1, 4) + sc*GJ*integ.t(x1,2) + d1*cos(theta_a)
    
    '''Row 8 v(x1) + sc*theta(x1) = d1cos(theta_a)'''
    GJ = (1/(G*(JJ.J())))
    FRzz = (1/(E*moi.I_zz()))
    
    ##A[8,0] = (FRzz/6) + (sc*GJ) #Ry1
    A[8,7] = -(FRzz)*x1  #Cv1
    A[8,8] = -(FRzz) #Cv2
    A[8,11] = sc*GJ #Ctheta
    '''RHS(8)'''
    b[8,0] = FRzz*integ.w(x1, 4) - (sc*GJ)*integ.t(x1,2) + d1*cos(theta_a)
    
    
    
    # '''Row 9 v(x2) - sc*theta(x2) = 0'''
    # A[9,0] = -(1/6)*((x2 - x1)**3)*(-FRzz) + sc*(x2-x1)*(-GJ)*sc #Ry1
    # A[9,6] = -(cos(theta_a)/6)*((x2 - (x2 - xa/2))**3)*(-FRzz) - cos(theta_a)*(0.5*ha - sc)*(x2 - (x2 - xa/2))*(-GJ)*sc + (sin(theta_a)*(ha/2)*(x2 - (x2 - xa/2))*(-GJ*sc))
    # A[9,7] = x2   #Cv1
    # A[9,8] = 1    #Cv2
    # A[9,11] = -sc #Ctheta
    # '''RHS(9)'''
    # b[9,0] = FRzz*integ.w(x2,4) + sc*GJ*(-P_y*(ha/2 - sc)*(x2 - (x2 - xa/2))) + sc*GJ*integ.t(x2,2)
    
    '''Row 9 v(x2) + sc*theta(x2) = 0'''
    A[9,0] = (FRzz/6)*((x2 - x1)**3) + (sc*GJ)*(x2-x1)*sc #Ry1
    ##A[9,2] = (FRzz/6) + (sc*GJ)*sc  #Ry2
    A[9,6] = (FRzz/6)*sin(theta_a)*((x2 - (x2 - xa/2))**3) + (sc*GJ)*(cos(theta_a)*(ha/2)*(x2 - (x2 - xa/2)) - (sin(theta_a)*(ha/2 - sc)*(x2 - (x2 - xa/2))))
    A[9,7] = -(FRzz)*x2   #Cv1
    A[9,8] = -(FRzz)    #Cv2
    A[9,11] = (sc*GJ) #Ctheta
    '''RHS(9)'''
    b[9,0] = FRzz*integ.w(x2,4) - (sc*GJ)*integ.t(x2,2)
    
    
    
    # '''Row 10 v(x3) - theta(x3)*sc = d3cos(theta_a)'''
    # A[10,0] = -(1/6)*((x3 - x1)**3)*(-FRzz) + sc*(x3 - x1)*(-GJ)*sc #Ry1
    # A[10,2] = -(1/6)*((x3 - x2)**3)*(-FRzz) + sc*(x3 - x2)*(-GJ)*sc #Ry2
    # A[10,6] = -(1/6)*sin(theta_a)*((x3 - (x2 - xa/2))**3)*(-FRzz) - sin(theta_a)*(ha/2 - sc)*(x3 - (x2 - xa/2))*(-GJ)*sc + cos(theta_a)*(ha/2)*(x3 - (x2 - xa/2))*(-GJ)*sc
    # A[10,7] = x3
    # A[10,9] = 1
    # A[10,11] = -1
    # '''RHS(10)'''
    # b[10, 0] = FRzz*(P_y/6)*((x3 - (x2 + xa/2))**3) + FRzz*integ.w(x3,4) + GJ*sc*((P_z*(ha/2)*(x3 - (x2 + xa/2)))- (P_y*(ha/2 - sc)*(x3 - (x2 - xa/2))) + (integ.t(x3,2))) + d3*cos(theta_a)
     
    '''Row 10 v(x3) + theta(x3)*sc = d3cos(theta_a)'''
    A[10,0] = (FRzz/6)*((x3 - x1)**3) + (sc*GJ)*(sc*(x3 - x1))  #Ry1
    A[10,2] = (FRzz/6)*((x3 - x2)**3) + (sc*GJ)*(sc*(x3 - x2))  #Ry2
    ##A[10,4] = (FRzz/6) + (sc*GJ)  #Ry3
    A[10,6] = (FRzz/6)*sin(theta_a)*((x3 - (x2 - xa/2))**3) + (sc*GJ)*(cos(theta_a)*(ha/2)*(x3 - (x2 - xa/2)) - sin(theta_a)*(ha/2 - sc)*(x3 - (x2 - xa/2))) 
    A[10,7] = -(FRzz)*x3 #Cv1
    A[10,8] = -(FRzz) #Cv2
    A[10,11] = (sc*GJ) #Ctheta
    '''RHS(10)'''
    b[10, 0] = (FRzz/6)*P_y*((x3 - (x2 + xa/2))**3) + FRzz*integ.w(x3,4) + (GJ*sc)*((P_z*(ha/2)*(x3 - (x2 + xa/2))) - (P_y*(ha/2 - sc)*(x3 - (x2 + xa/2))) - (integ.t(x3,2))) + d3*cos(theta_a)
    
    
    
    # '''Row 11 w(x2 + xa/2)cos + v(x2 + xa/2)sin'''
    # A[11,0] = -(1/6)*(((x2 - xa/2) - x1)**3)*(-FRzz)*sin(theta_a)
    # A[11,1] = -(1/6)*(((x2 - xa/2) - x1)**3)*(-FRyy)*cos(theta_a)
    # A[11,7] = (x2 - xa/2)*sin(theta_a)
    # A[11,8] = sin(theta_a)
    # A[11,9] = (x2 - xa/2)*cos(theta_a)
    # A[11,10] = cos(theta_a)
    # '''RHS(11)'''
    # b[11,0] = FRzz*integ.w((x2 - xa/2), 4)
    

    '''Row 11 (v(xaI) - theta(xaI)*(ha/2 - sc))*sin(theta_a) + (w(xaI) + theta(xaI)*(ha/2))*cos(theta_a) = 0'''    
    x = x2 - xa/2
    
    A[11, 0] = ((FRzz/6)*((x - x1)**3) - (GJ*sc)*(ha/2 - sc)*(x - x1))*sin(theta_a) + ((GJ*sc)*(ha/2)*((x - x1)**3))*cos(theta_a)
    A[11, 1] = (FRyy/6)*((x - x1)**3)*cos(theta_a)
    ##A[11, 6] = ((FRzz/6)*sin(theta_a) - GJ*(cos(theta_a)*(ha/2) - (ha/2 - sc)*sin(theta_a)))*sin(theta_a) + (FRyy*cos(theta_a) + GJ*((ha/2)*cos(theta_a) - (ha/2 - sc)*sin(theta_a)))*cos(theta_a)
    A[11, 7] = (-FRzz*x)*sin(theta_a) 
    A[11, 8] = (-FRzz)*sin(theta_a)
    A[11, 9] = (-FRyy*x)*cos(theta_a)
    A[11,10] = (-FRyy)*cos(theta_a)
    A[11,11] = (-GJ)*(ha/2 - sc)*sin(theta_a) + GJ*(ha/2)*cos(theta_a)   
    '''RHS(11)'''
    b[11, 0] = (FRzz*integ.w(x, 4) + GJ*(ha/2 - sc)*integ.t(x, 2))*sin(theta_a) - (GJ*(ha/2)*integ.t(x, 2))*cos(theta_a)
    
    
    
    
    xx = np.linalg.solve(A, b)
    
    Fy = -integ.w(la,1) - P_y
    Fz = -P_z
    
    return (A, xx, Fy, Fz)




# AA, xx, FY, FZ= getunknowns()

# print ("Ry", xx[0] + xx[2] + xx[4] + xx[6]*0.43837114)
# print ("Fy", FY)
# print ("Rz", xx[1] + xx[3] + xx[5] + xx[6]*0.89879405)
# print ("FZ", FZ)









