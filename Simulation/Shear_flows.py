import numpy as np
import Num_integration as Nint
import Coordinates
from numpy import nextafter as nf
import GlobalConstants as G
import MoI
import Centroid as Ce

#These functions return the shear flows that a unit vertical shear, horizontal shear, and pure torque casue, respectively
#you can then use "stress_plots.py" to plot them
def SF_vertical(Ca = G.Ca, h=G.ha, tsk=G.tsk, tsp=G.tsp, tst=G.tst, hst=G.hst, wst=G.wst, Nst=G.nst, Izz=MoI.I_zz(), B737 = 0):  ### Chord length, aileron height, coordinates of stringers, skin thickness,
    ### spar thickness, stringer #, stringer area, I about z-axis (should be previously calculated)
    Coord = Coordinates.Coord_out(Ca, h, Nst)
    Ast = hst*tst + wst*tst
    r = h/2
    Per = np.pi * r + 2 * np.sqrt((Ca - r) ** 2 + r ** 2)  # perimeter
    Dis = Per / Nst  # distance between stringers
    theta = Dis / r  # Angle between first stringer and the adjacent ones, in radians
    lsk = np.sqrt((Ca - r) ** 2 + r ** 2)  # Length of diagonal sections of skin of aileron
    smalldist = Dis - (np.pi/2 - theta) * r #Small distance between junction of spar and skin and the next adjacent boom

    #First calculating the base shear flows in each section due to a unit load applied at the S.C.:
    #Section 1
    def qb_LE_top(th):
        if th < theta:
            qb = (-1/Izz)*tsk*r**2*Nint.Simpson38_int(np.sin, 0, th)
        else:
            qb = (-1/Izz)*(tsk*r**2*Nint.Simpson38_int(np.sin, 0, th) + Ast*Coord[1, 1])
        return qb
    #Section 2
    def qb_sp_top(s):
        qb = (-1 / Izz) * tsp * Nint.Simpson_int(lambda y: y, 0, s)
        return qb
    #print(qb_sp_top(r))

    #Section 3
    if B737 == 0:
        def qb_sk_top(s):
            if s < smalldist:
                qb = (-1 / Izz) * Nint.Simpson_int(lambda y: tsk*(r - (r/(lsk))*y), 0, s) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: tsk*(r - (r/(lsk))*y), 0, s) + Ast * Coord[2, 1]) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 2*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: tsk*(r - (r/(lsk))*y), 0, s) + Ast * np.sum(Coord[2:4, 1])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 3*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: tsk*(r - (r/(lsk))*y), 0, s) + Ast * np.sum(Coord[2:5, 1])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 4*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: tsk*(r - (r/(lsk))*y), 0, s) + Ast * np.sum(Coord[2:6, 1])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            else:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: tsk*(r - (r/(lsk))*y), 0, s) + Ast * np.sum(Coord[2:7, 1])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            return qb

    else:
        def qb_sk_top(s):
            if s < smalldist:
                qb = (-1 / Izz) * Nint.Simpson_int(lambda y: tsk*(r - (r/(lsk))*y), 0, s) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: tsk*(r - (r/(lsk))*y), 0, s) + Ast * Coord[2, 1]) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 2*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: tsk*(r - (r/(lsk))*y), 0, s) + Ast * np.sum(Coord[2:4, 1])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 3*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: tsk*(r - (r/(lsk))*y), 0, s) + Ast * np.sum(Coord[2:5, 1])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 4*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: tsk*(r - (r/(lsk))*y), 0, s) + Ast * np.sum(Coord[2:6, 1])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 5*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: tsk*(r - (r/(lsk))*y), 0, s) + Ast * np.sum(Coord[2:7, 1])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            else:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: tsk*(r - (r/(lsk))*y), 0, s) + Ast * np.sum(Coord[2:8, 1])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            return qb

    #Section 4
    if B737 ==0:
        def qb_sk_bot(s):
            if s < Dis/2:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: -(r/(lsk))*y*tsk, 0, s)) + qb_sk_top(lsk)
            elif s < Dis/2 + Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: -(r/(lsk))*y*tsk, 0, s) + Ast * Coord[7, 1]) + qb_sk_top(lsk)
            elif s < Dis/2 + 2*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: -(r/(lsk))*y*tsk, 0, s) + Ast * np.sum(Coord[7:9, 1])) + qb_sk_top(lsk)
            elif s < Dis/2 + 3*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: -(r/(lsk))*y*tsk, 0, s) + Ast * np.sum(Coord[7:10, 1])) + qb_sk_top(lsk)
            elif s < Dis/2 + 4*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: -(r/(lsk))*y*tsk, 0, s) + Ast * np.sum(Coord[7:11, 1])) + qb_sk_top(lsk)
            else:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: -(r/(lsk))*y*tsk, 0, s) + Ast * np.sum(Coord[7:12, 1])) + qb_sk_top(lsk)
            return qb
    else:
        def qb_sk_bot(s):
            if s < Dis/2:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: -(r/(lsk))*y*tsk, 0, s)) + qb_sk_top(lsk)
            elif s < Dis/2 + Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: -(r/(lsk))*y*tsk, 0, s) + Ast * Coord[8, 1]) + qb_sk_top(lsk)
            elif s < Dis/2 + 2*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: -(r/(lsk))*y*tsk, 0, s) + Ast * np.sum(Coord[8:10, 1])) + qb_sk_top(lsk)
            elif s < Dis/2 + 3*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: -(r/(lsk))*y*tsk, 0, s) + Ast * np.sum(Coord[8:11, 1])) + qb_sk_top(lsk)
            elif s < Dis/2 + 4*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: -(r/(lsk))*y*tsk, 0, s) + Ast * np.sum(Coord[8:12, 1])) + qb_sk_top(lsk)
            elif s < Dis/2 + 5*Dis:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: -(r/(lsk))*y*tsk, 0, s) + Ast * np.sum(Coord[8:13, 1])) + qb_sk_top(lsk)
            else:
                qb = (-1 / Izz) * (Nint.Simpson_int(lambda y: -(r/(lsk))*y*tsk, 0, s) + Ast * np.sum(Coord[8:14, 1])) + qb_sk_top(lsk)
            return qb

    #Section 5
    def qb_sp_bot(s):
        qb = (-1 / Izz) * tsp * Nint.Simpson_int(lambda y: y, 0, s)
        return qb
    #print(qb_sp_bot(-r))

    #Section 6
    def qb_LE_bot(th):
        if th < -theta:
            qb = (-1/Izz)*tsk*r**2*Nint.Simpson38_int(np.sin, -np.pi/2, th) + qb_sk_bot(lsk) - qb_sp_bot(-r)
        else:
            qb = (-1/Izz)*(tsk*r**2*Nint.Simpson38_int(np.sin, -np.pi/2, th) + Ast*Coord[(Coord.shape[0]-1), 1]) + qb_sk_bot(lsk) - qb_sp_bot(-r)
        return qb

    #Zero twist condition to calculate redundant shear flows, done for each "cell" separately
    #Calculations for the right cell, cell I:

    #Contribution from leading edge, section 1 and section 6:
    int_qb1 = r*(Nint.Simpson38_int(qb_LE_top, 0, np.nextafter(theta, 0)) + Nint.Simpson38_int(qb_LE_top, np.nextafter(theta, theta+1), np.pi/2) )
    int_qb6 = r*( Nint.Simpson38_int(qb_LE_bot, -np.pi/2, np.nextafter(-theta, -theta-1)) + Nint.Simpson38_int(qb_LE_bot, np.nextafter(-theta, -theta+1), 0)  )
    #Contribution from spar, section 2 and 5:
    int_qb2I = -Nint.Simpson_int(qb_sp_top, r, 0)
    int_qb5I = -Nint.Simpson_int(qb_sp_bot, 0, -r)

    int_cellI = int_qb1/tsk + int_qb6/tsk + int_qb2I/tsp + int_qb5I/tsp

    #Calculations for the left cell, cellII:
    #Contribution from top and bottom skin sections, section 3 and section 4:
    sd = smalldist
    int_qb3 = Nint.Simpson_int(qb_sk_top, 0, np.nextafter(sd, 0)) + Nint.Simpson_int(qb_sk_top, nf(sd, 1), nf(sd + Dis, 0)) \
              + Nint.Simpson_int(qb_sk_top, nf(sd + Dis, 1), nf(sd + 2*Dis, 0)) + Nint.Simpson_int(qb_sk_top, nf(sd + 2*Dis, 1), nf(sd + 3*Dis, 0)) \
              + Nint.Simpson_int(qb_sk_top, nf(sd + 3*Dis, 1), nf(sd + 4*Dis, 0)) + Nint.Simpson_int(qb_sk_top, nf(sd + 4*Dis, 1), lsk)

    int_qb4 = Nint.Simpson_int(qb_sk_bot, 0, np.nextafter(Dis/2, 0)) + Nint.Simpson_int(qb_sk_bot, nf(Dis/2, 1), nf(1.5*Dis, 0)) \
              + Nint.Simpson_int(qb_sk_bot, nf(1.5*Dis, 1), nf(2.5*Dis, 0)) + Nint.Simpson_int(qb_sk_bot, nf(2.5*Dis, 1), nf(3.5*Dis, 0)) \
              + Nint.Simpson_int(qb_sk_bot, nf(3.5*Dis, 1), nf(4.5*Dis, 0)) + Nint.Simpson_int(qb_sk_bot, nf(4.5*Dis, 1), lsk)

    #Contribution from spar, section 2 and 5:
    int_qb5II = Nint.Simpson_int(qb_sp_bot, -r, 0)
    int_qb2II = Nint.Simpson_int(qb_sp_top, 0, r)

    int_cellII = int_qb3/tsk + int_qb4/tsk + int_qb5II/tsp + int_qb2II/tsp

    x11 = (((r * np.pi) / tsk) + h / tsp)
    x21 = ((-h / tsp))
    x12 = ((-h / tsp))
    x22 = (((2 * lsk) / tsk) + h / tsp)
    c1 = -int_cellI
    c2 = -int_cellII

    a = np.array([[x11, x21], [x12, x22]])
    b = np.array([c1, c2])
    sol = np.linalg.solve(a, b)
    #print(sol)

    #New functions for total shear flow (base plus calculated redundants)
    def q1_LE_top(th):
        return qb_LE_top(th) + sol[0]
    def q2_sp_top(s):
        return qb_sp_top(s) - sol[0] + sol[1]
    def q3_sk_top(s):
        return qb_sk_top(s) + sol[1]
    def q4_sk_bot(s):
        return qb_sk_bot(s) + sol[1]
    def q5_sp_bot(s):
        return qb_sp_bot(s) - sol[0] + sol[1]
    def q6_LE_bot(th):
        return qb_LE_bot(th) + sol[0]

    return (q1_LE_top, q2_sp_top, q3_sk_top, q4_sk_bot, q5_sp_bot, q6_LE_bot)

def SF_horizontal(Ca = G.Ca, h=G.ha, tsk=G.tsk, tsp=G.tsp, tst=G.tst, hst=G.hst, wst=G.wst, Nst=G.nst, Iyy=MoI.I_yy(), B737 = 0):
    Coord = Coordinates.Coord_out(Ca, h, Nst)
    zsc = -Ce.Centroid_z() #defining centroid distance as a positive distance from spar, to the right
    Coord[:, 0] = Coord[:, 0] + zsc #transforming coordinates of booms to a centroidal coordinate system
    Ast = hst*tst + wst*tst
    r = h/2
    Per = np.pi * r + 2 * np.sqrt((Ca - r) ** 2 + r ** 2)  # perimeter
    Dis = Per / Nst  # distance between stringers
    theta = Dis / r  # Angle between first stringer and the adjacent ones, in radians
    lsk = np.sqrt((Ca - r) ** 2 + r ** 2)  # Length of diagonal sections of skin of aileron
    smalldist = Dis - (np.pi/2 - theta) * r

    #Section 1
    def qb_LE_top(th):
        if th < theta:
            qb = (-1/Iyy)*(tsk*r*Nint.Simpson38_int(lambda th1: r*np.cos(th1) + zsc, 0, th) + 0.5*Ast*Coord[0, 0]) #Note the half area term due to cut at boom
        else:
            qb = (-1/Iyy)*(tsk*r*Nint.Simpson38_int(lambda th1: r*np.cos(th1) + zsc, 0, th) + 0.5*Ast*Coord[0, 0] + Ast*Coord[1, 0])
        return qb

    #Section 2
    def qb_sp_top(s):
        qb = (-1 / Iyy) * tsp * zsc*s
        return qb
    #print(qb_sp_top(r))

    #Section 3
    if B737 == 0:
        def qb_sk_top(s):
            if s < smalldist:
                qb = (-1 / Iyy) * Nint.Simpson_int(lambda s3: tsk*(zsc - ((Ca-r)/(lsk))*s3), 0, s) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s3: tsk*(zsc - ((Ca-r)/(lsk))*s3), 0, s) + Ast * Coord[2, 0]) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 2*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s3: tsk*(zsc - ((Ca-r)/(lsk))*s3), 0, s) + Ast * np.sum(Coord[2:4, 0])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 3*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s3: tsk*(zsc - ((Ca-r)/(lsk))*s3), 0, s) + Ast * np.sum(Coord[2:5, 0])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 4*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s3: tsk*(zsc - ((Ca-r)/(lsk))*s3), 0, s) + Ast * np.sum(Coord[2:6, 0])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            else:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s3: tsk*(zsc - ((Ca-r)/(lsk))*s3), 0, s) + Ast * np.sum(Coord[2:7, 0])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            return qb

    else:
        def qb_sk_top(s):
            if s < smalldist:
                qb = (-1 / Iyy) * Nint.Simpson_int(lambda s3: tsk*(zsc - ((Ca-r)/(lsk))*s3), 0, s) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s3: tsk*(zsc - ((Ca-r)/(lsk))*s3), 0, s) + Ast * Coord[2, 0]) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 2*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s3: tsk*(zsc - ((Ca-r)/(lsk))*s3), 0, s) + Ast * np.sum(Coord[2:4, 0])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 3*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s3: tsk*(zsc - ((Ca-r)/(lsk))*s3), 0, s) + Ast * np.sum(Coord[2:5, 0])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 4*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s3: tsk*(zsc - ((Ca-r)/(lsk))*s3), 0, s) + Ast * np.sum(Coord[2:6, 0])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            elif s < smalldist + 5*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s3: tsk*(zsc - ((Ca-r)/(lsk))*s3), 0, s) + Ast * np.sum(Coord[2:7, 0])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            else:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s3: tsk*(zsc - ((Ca-r)/(lsk))*s3), 0, s) + Ast * np.sum(Coord[2:8, 0])) + qb_sp_top(r) + qb_LE_top(np.pi/2)
            return qb

    #Section 4
    if B737 ==0:
        def qb_sk_bot(s):
            if s < Dis/2:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s4: tsk*(zsc - (Ca-r) + s4*(Ca-r)/lsk), 0, s)) + qb_sk_top(lsk)
            elif s < Dis/2 + Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s4: tsk*(zsc - (Ca-r) + s4*(Ca-r)/lsk), 0, s) + Ast * Coord[7, 0]) + qb_sk_top(lsk)
            elif s < Dis/2 + 2*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s4: tsk*(zsc - (Ca-r) + s4*(Ca-r)/lsk), 0, s) + Ast * np.sum(Coord[7:9, 0])) + qb_sk_top(lsk)
            elif s < Dis/2 + 3*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s4: tsk*(zsc - (Ca-r) + s4*(Ca-r)/lsk), 0, s) + Ast * np.sum(Coord[7:10, 0])) + qb_sk_top(lsk)
            elif s < Dis/2 + 4*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s4: tsk*(zsc - (Ca-r) + s4*(Ca-r)/lsk), 0, s) + Ast * np.sum(Coord[7:11, 0])) + qb_sk_top(lsk)
            else:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s4: tsk*(zsc - (Ca-r) + s4*(Ca-r)/lsk), 0, s) + Ast * np.sum(Coord[7:12, 0])) + qb_sk_top(lsk)
            return qb
    else:
        def qb_sk_bot(s):
            if s < Dis/2:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s4: tsk*(zsc - (Ca-r) + s4*(Ca-r)/lsk), 0, s)) + qb_sk_top(lsk)
            elif s < Dis/2 + Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s4: tsk*(zsc - (Ca-r) + s4*(Ca-r)/lsk), 0, s) + Ast * Coord[8, 0]) + qb_sk_top(lsk)
            elif s < Dis/2 + 2*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s4: tsk*(zsc - (Ca-r) + s4*(Ca-r)/lsk), 0, s) + Ast * np.sum(Coord[8:10, 0])) + qb_sk_top(lsk)
            elif s < Dis/2 + 3*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s4: tsk*(zsc - (Ca-r) + s4*(Ca-r)/lsk), 0, s) + Ast * np.sum(Coord[8:11, 0])) + qb_sk_top(lsk)
            elif s < Dis/2 + 4*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s4: tsk*(zsc - (Ca-r) + s4*(Ca-r)/lsk), 0, s) + Ast * np.sum(Coord[8:12, 0])) + qb_sk_top(lsk)
            elif s < Dis/2 + 5*Dis:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s4: tsk*(zsc - (Ca-r) + s4*(Ca-r)/lsk), 0, s) + Ast * np.sum(Coord[8:13, 0])) + qb_sk_top(lsk)
            else:
                qb = (-1 / Iyy) * (Nint.Simpson_int(lambda s4: tsk*(zsc - (Ca-r) + s4*(Ca-r)/lsk), 0, s) + Ast * np.sum(Coord[8:14, 0])) + qb_sk_top(lsk)
            return qb

    #Section 5
    def qb_sp_bot(s):
        qb = (-1 / Iyy) * tsp * zsc * s
        return qb
    #print(qb_sp_bot(-r))

    #Section 6
    def qb_LE_bot(th):
        if th < -theta:
            qb = (-1/Iyy)*tsk*r*Nint.Simpson38_int(lambda th1: r*np.cos(th1) + zsc, -np.pi/2, th) + qb_sk_bot(lsk) - qb_sp_bot(-r)
        elif th == 0:
            qb = (-1/Iyy)*(tsk*r*Nint.Simpson38_int(lambda th1: r*np.cos(th1) + zsc, -np.pi/2, th) + 0.5*Ast*Coord[0,0] +
                           Ast*Coord[(Coord.shape[0]-1), 0]) + qb_sk_bot(lsk) - qb_sp_bot(r)
        else:
            qb = (-1/Iyy)*(tsk*r*Nint.Simpson38_int(lambda th1: r*np.cos(th1) + zsc, -np.pi/2, th) + Ast*Coord[(Coord.shape[0]-1), 0]) + qb_sk_bot(lsk) - qb_sp_bot(-r)
        return qb
#Redundant shear flows are zero due to cuts at symmetry axis

    return (qb_LE_top, qb_sp_top, qb_sk_top, qb_sk_bot, qb_sp_bot, qb_LE_bot)

def SF_torque (Ca = G.Ca, ha = G.ha, tsk = G.tsk, tsp = G.tsp):
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
    return(q01, q02)
