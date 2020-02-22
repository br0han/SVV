import GlobalConstants as GC
import numpy as np
from matplotlib import pyplot as plt
from cubicspline import cubicspline
from aerodynamicloads import DistForceTorqueMatrix
"""
data = np.genfromtxt('aerodynamicloadcrj700.dat', dtype=None, delimiter=',')
val = data[:,1]
ShearCenter = 0.2

n=len(val)
Ca = GC.Ca #0.484
la = GC.la #1.691
Coord = []
for i in range(1,n+1):
    theta0 = (i-1)/n * np.pi
    theta1 = (i)/n * np.pi
    z = -1/2*(Ca/2*(1-np.cos(theta0))+Ca/2*(1-np.cos(theta1)))
    Coord.append(z)
Coord = np.array(Coord)
"""
MatrixF , MatrixT , MatrixMF, MatrixMT = DistForceTorqueMatrix() #data,Coord,ShearCenter

<<<<<<< HEAD
=======
MatrixF , MatrixT , MatrixMF, MatrixMT = DistForceTorqueMatrix()
>>>>>>> 9c13501ddf095838e1e6df9f27ff586170a31694

n = len(MatrixF)+1
Ca = GC.Ca #0.484
la = GC.la #1.691
Span = []
for i in range(1,n+1):
    theta0 = (i-1)/n * np.pi
    theta1 = (i)/n * np.pi
    z = 1/2*(la/2*(1-np.cos(theta0))+la/2*(1-np.cos(theta1)))
    Span.append(z)
Span = np.array(Span)



XX = list(np.linspace(Span[0],Span[-1],10000))
FF = []
for i in np.linspace(Span[0],Span[-1],10000):
    for j in range(0,len(MatrixF)):
        if Span[j+1] >=  i >= Span[j]:
            
            a = MatrixF[j,0]
            b = MatrixF[j,1]
            c = MatrixF[j,2]
            d = MatrixF[j,3]

            y = a * ( i - Span[j])**3 + b*(i-Span[j])**2 + c*(i-Span[j]) + d
            

            FF.append(y)

       
            
XX = list(np.linspace(Span[0],Span[-1],10000))
TT = []
for i in np.linspace(Span[0],Span[-1],10000):
    for j in range(0,len(MatrixT)):
        if Span[j+1] >=  i >= Span[j]:
            
            a = MatrixT[j,0]
            b = MatrixT[j,1]
            c = MatrixT[j,2]
            d = MatrixT[j,3]

            y = a * ( i - Span[j])**3 + b*(i-Span[j])**2 + c*(i-Span[j]) + d
            

            TT.append(y)


            
plt.subplot(1,2,1)
plt.plot(XX,FF)
plt.legend('Force')

'''plt.subplot(1,2,2)
plt.plot(XX,TT)
plt.legend('Torque')
'''
plt.grid(1)

plt.show()

