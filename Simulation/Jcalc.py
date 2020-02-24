import GlobalConstants as Gc
import numpy as np
import Coordinates

def J (Ca = Gc.Ca, ha = Gc.ha, tsk = Gc.tsk, tsp = Gc.tsp):
    r = ha/2
    lsk = np.sqrt((Ca - r) ** 2 + r ** 2)
    A1 = np.pi*(r**2)*0.5
    A2 = r*(Ca-r)

    #Coefficients for equations
    x11 = 2*A1
    x12 = 2*A2
    x21 = (1/(2*A1))*( (np.pi*r)/tsk + ha/tsp) + (ha/tsp)*(1/(2*A2))
    x22 = -(1/(2*A1))*(ha/tsp) - (1/(2*A2))*( (2*lsk)/tsk + ha/tsp)

    #The constants the equations have to be equal to, first one is 1 due to the assumed unit torque
    c1 = 1
    c2 = 0

    a = np.array([[x11, x12], [x21, x22]])
    b = np.array([c1, c2])
    q01, q02 = np.linalg.solve(a, b)

    J = 1 / ((1/(2*A1))*( (q01*np.pi*r)/tsk + (q01 - q02)*(ha/tsp) ))
    return(J)


