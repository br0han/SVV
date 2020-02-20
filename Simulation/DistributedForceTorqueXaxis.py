import GlobalConstants as GC
import numpy as np
from matplotlib import pyplot as plt
from cubicspline import cubicspline
from aerodynamicloads import CoordForceTorque

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

DT , DF = CoordForceTorque(data,Coord,ShearCenter)

n = len(DT)
Ca = GC.Ca #0.484
la = GC.la #1.691
Span = []
for i in range(1,n+1):
    theta0 = (i-1)/n * np.pi
    theta1 = (i)/n * np.pi
    z = -1/2*(la/2*(1-np.cos(theta0))+la/2*(1-np.cos(theta1)))
    Span.append(z)
Span = np.array(Span)

MatrixF = cubicspline(DT , Span)
MatrixT = cubicspline(DF , Span)

XX = list(np.linspace(Span[0],Span[-1],10000))
YY = []
for i in np.linspace(Span[0],Span[-1],10000):
    for j in range(0,39):
        if Span[j+1] <=  i <= Span[j]:
            
            a = MatrixT[j,0]
            b = MatrixT[j,1]
            c = MatrixT[j,2]
            d = MatrixT[j,3]

            #a = cubicspline(val,Coord)[j,0]
            #b = cubicspline(val,Coord)[j,1]
            #c = cubicspline(val,Coord)[j,2]
            #d = cubicspline(val,Coord)[j,3]
            

            y = a * ( i - Span[j])**3 + b*(i-Span[j])**2 + c*(i-Span[j]) + d
            

            YY.append(y)


            


plt.plot(XX,YY)
plt.show()

