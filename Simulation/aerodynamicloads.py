import GlobalConstants as GC
import numpy as np
from matplotlib import pyplot as plt
from cubicspline import cubicspline


data = np.genfromtxt('aerodynamicloadcrj700.dat', dtype=None, delimiter=',')
val = data[:,1]

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


CoefM = cubicspline(val,Coord)


