import Coordinates as f
import GlobalConstants as g
from math import *

def Centroid_z():

    C_a = g.Ca     #[m]
    h = g.ha       #[m]
    t_sk = g.tsk   #[m]
    t_sp = g.tsp   #[m]
    h_st = g.hst    #[m]
    w_st = g.wst    #[m]
    t_st = g.tst

    R = h/2

    Coord_z = (f.Coord_out()[:,0])

    #[half circle centroid]
    z_1 = (2*R)/pi
    a_1 = pi*R*t_sk

    #[triangle]
    z_2 = -0.5*(C_a-R)
    a_2 = sqrt((C_a-R)**2+(R**2))*t_sk

    #[spar]
    z_3 = 0
    a_3 = t_sp*h

    #[stiddeners]
    z_4 = sum(Coord_z)
    a_4 = (h_st + w_st)*t_st

    #[centroid]
    z = ((z_1*a_1)+(2*(z_2*a_2))+(z_3*a_3)+(z_4*a_4))/(a_1+(2*a_2)+a_3+(13*a_4))

    return(z)

