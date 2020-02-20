import GlobalConstants as GC
import numpy as np
from numpy.linalg import inv
from matplotlib import pyplot as plt

#INPUTS: val   = values of the function i.e. data points, inputed as a numpy array
#        Coord = the coordinates of the data points, inputed as a numpy array

def cubicspline(val, Coord):

    #--Coefficients of M matrix
    n=len(val)
    matrix = np.zeros((n-2,n-2))

    for j in range(0,n-2):
        Cof1 = (Coord[j+1]-Coord[j])/6
        Cof2 = ((Coord[j+1]-Coord[j])+(Coord[j+2]-Coord[j+1]))/3
        Cof3 = (Coord[j+2]-Coord[j+1])/6

        if j == 0:
            matrix[j,0] = Cof2
            matrix[j,1] = Cof3
        
        elif j == n-3:
            matrix[j,n-4] = Cof1
            matrix[j,n-3] = Cof2
        else:
            matrix[j,j-1]  = Cof1
            matrix[j,j]    = Cof2
            matrix[j,j+1]  = Cof3

    #--RHS column            
    Col = np.zeros(n-2)

    for j in range (0,n-2):
        Col[j] = (val[j+2]-val[j+1])/(Coord[j+2]-Coord[j+1]) - (val[j+1]-val[j])/(Coord[j+1]-Coord[j])

    Col = Col.reshape((n-2,1))


    #--Calculating M values
    M = np.matmul(inv(matrix),Col)

    #--Setting Boundary Condition (in this case, natural spline)
    M0 = np.array([0])
    Mn = np.array([0])

    M = np.vstack((M0,M))
    M = np.vstack((M,Mn))

    #--Calculating Spline Coefficients a,b,c,d
    CoefM = np.zeros((n-1,4))

    for j in range(0,n-1):
        M_i   = M[j,0]
        M_i1  = M[j+1,0]
        h_i   = Coord[j+1]-Coord[j]
        f_i   = val[j]
        f_i1 = val[j+1]
    
    
        CoefM[j,0] = (M_i1 - M_i)/(6 * h_i)
        CoefM[j,1] =  M_i / 2
        CoefM[j,2] = (f_i1 - f_i)/h_i - h_i / 3 * M_i - h_i/6* M_i1
        CoefM[j,3] = f_i
    
    return CoefM   #4 colomuns: a, b, c, d respectively. Rows are equations for different spline


















