import Shear_flows as SF
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import math
import Centroid as Ce
import GlobalConstants as G
import MoI

r = G.ha/2
zsc = -Ce.Centroid_z()
Ca = G.Ca
lsk = np.sqrt((Ca - r) ** 2 + r ** 2)

#The get data functions use the functions from "Shear_flows.py" to arrange the data in a proper way to input it into the
#plotting function. Each of them outputs a n by 3 array, first column is shear flow, second z-coordinate, last is y.
#The input scales the flows properly with the given load.

def getSFDataTot(Sy = 1, Sz = 1, T = 1, FlowOrStress = 1): #Adds all the shear flows properly. The last parameter determines
                                                           #whether shear FLOW or shear STRESS is outputted, outputs flows by default,
                                                           #change the 'FlowOrStress" to not one for stress outputting
    tsk = G.tsk
    tsp = G.tsp
    if FlowOrStress == 1:
        tsk = 1
        tsp = 1
    else:
        tsk = tsk
        tsp = tsp

    q1h, q2h, q3h, q4h, q5h, q6h = SF.SF_horizontal()
    q1v, q2v, q3v, q4v, q5v, q6v = SF.SF_vertical()
    q01, q02 = SF.SF_torque()
    q01, q02 = -q01, -q02 #Clockwise shear flow is positive, but counterlockwise torque is positive

    data = np.zeros((1, 3))

    for th in np.arange(0, np.pi/2, np.pi/200):
        if th == 0:
            data[0, 0] = (Sy*q1v(th) + Sz*q1h(th) + T*q01)/tsk
            data[0, 1] = r*np.cos(th) + zsc
            data[0, 2] = r*np.sin(th)
        else:
            data = np.append(data, [[(Sy*q1v(th) + Sz*q1h(th) + T*q01)/tsk, r*np.cos(th) + zsc, r*np.sin(th)]], axis=0)
    for s2 in np.arange(0, r, r/80):
        data = np.append(data, [[(Sy*q2v(s2) + Sz*q2h(s2) + T*(q02 - q01))/tsp, zsc, s2]], axis=0)
    for s3 in np.arange(0, lsk, lsk/400):
        data = np.append(data, [[(Sy*q3v(s3) + Sz*q3h(s3) + T*q02)/tsk, zsc - ((Ca-r)/(lsk))*s3, r - (r/(lsk))*s3]], axis=0)
    for s4 in np.arange(0, lsk, lsk/400):
        data = np.append(data, [[(Sy*q4v(s4) + Sz*q4h(s4) + T*q02)/tsk, zsc - (Ca-r) + s4*((Ca-r)/lsk), -(r/(lsk))*s4]], axis=0)
    for s5 in np.arange(0, -r, -r/80):
        data = np.append(data, [[(Sy*q5v(s5) + Sz*q5h(s5) + T*(q02 - q01))/tsp, zsc, s5]], axis=0)
    for th in np.arange(-np.pi/2+0.025, 0, (np.pi + 0.025)/200):
        data = np.append(data, [[(Sy*q6v(th) + Sz*q6h(th) + T*q01)/tsk, r*np.cos(th) + zsc, r*np.sin(th)]], axis=0)
    return data

def getNSDataTot(My = 1, Mz = 1, Izz = MoI.I_zz(), Iyy = MoI.I_yy()): #Adds all the normal stresses from the moments

    data = np.zeros((1, 3))

    for th in np.arange(0, np.pi/2, np.pi/200):
        z = r*np.cos(th) + zsc
        y = r*np.sin(th)
        if th == 0:
            data[0, 0] = (Mz*y)/Izz + (My*z)/Iyy
            data[0, 1] = z
            data[0, 2] = y
        else:
            data = np.append(data, [[(Mz*y)/Izz + (My*z)/Iyy, r*np.cos(th) + zsc, r*np.sin(th)]], axis=0)
    for s2 in np.arange(0, r, r/80):
        z = zsc
        y = s2
        data = np.append(data, [[(Mz*y)/Izz + (My*z)/Iyy, z, y]], axis=0)
    for s3 in np.arange(0, lsk, lsk/400):
        z = zsc - ((Ca-r)/(lsk))*s3
        y = r - (r/(lsk))*s3
        data = np.append(data, [[(Mz*y)/Izz + (My*z)/Iyy, z, y]], axis=0)
    for s4 in np.arange(0, lsk, lsk/400):
        z =  zsc - (Ca-r) + s4*((Ca-r)/lsk)
        y = -(r/(lsk))*s4
        data = np.append(data, [[(Mz*y)/Izz + (My*z)/Iyy, z, y]], axis=0)
    for s5 in np.arange(0, -r, -r/80):
        z = zsc
        y = s5
        data = np.append(data, [[(Mz*y)/Izz + (My*z)/Iyy, z, y]], axis=0)
    for th in np.arange(-np.pi/2+0.025, 0, (np.pi + 0.025)/200):
        z = r*np.cos(th) + zsc
        y = r*np.sin(th)
        data = np.append(data, [[(Mz*y)/Izz + (My*z)/Iyy, r*np.cos(th) + zsc, r*np.sin(th)]], axis=0)
    return data

def getVMData(Sy=1, Sz=1, My=1, Mz=1, T=1):
    dataNS = getNSDataTot(My, Mz)
    dataSF = getSFDataTot(Sy, Sz, T, FlowOrStress=0)
    dataVM = np.zeros(dataNS.shape)
    dataVM[:, 1:3] = dataNS[:, 1:3]
    dataVM[:, 0] = np.sqrt((dataNS[:, 0])**2 + 3*(dataSF[:, 0])**2 +0.000000000001)
    return dataVM

def plot ( data, fname, label, x):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    plt.xlim(0.225, -0.325)
    plt.ylim(-0.125, 0.125)
    plt.xlabel("z [m]")
    plt.ylabel("y [m]")
    plt.scatter(data[:, 1], data[:, 2], c=data[:, 0], s=3)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    plt.colorbar(cax=cax, label=label)
    plt.set_cmap("jet")
    plt.show()
plot(getSFDataTot(1, 0, 0), "hello", 'label', 2)