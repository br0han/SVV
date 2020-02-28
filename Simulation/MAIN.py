# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 22:37:06 2020

@author: Burhan

MAIN
"""

import time

t0 = time.time()


import numpy as np
import GlobalConstants as g
import InternalMomentsAndDeflections as graph
import matplotlib.pyplot as plt

N = 10000 #number of points to sample

x = np.linspace(0, g.la, N)

Mz = np.zeros_like(x)
My = np.zeros_like(x)
Sy = np.zeros_like(x)
Sz = np.zeros_like(x)
T = np.zeros_like(x)
v = np.zeros_like(x)
w = np.zeros_like(x)
theta = np.zeros_like(x)



for i in range(N):
    
    My[i] = graph.My(x[i])
    Mz[i] = graph.Mz(x[i])
    Sy[i] = graph.Sy(x[i])
    Sz[i] = graph.Sz(x[i])
    T[i] = graph.T(x[i])
    v[i] = graph.v(x[i])
    w[i] = graph.w(x[i])
    theta[i] = graph.theta(x[i])
    
print ("That took", time.time() - t0, "second(s)")


'''#####################################################################################'''
'bo'
fig1, ax1 = plt.subplots(1, figsize = (8, 8))
ax1.plot(x, My, linewidth=4);
ax1.set_title('Moment in y direction');
ax1.grid(True);
ax1.set_xlabel("x [m]");
ax1.set_ylabel("My [kNm]");
plt.tight_layout()

fig2, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, Mz, linewidth=4);
ax.set_title('Moment in z direction');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("Mz [kNm]");

fig3, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, Sy, linewidth=4);
ax.set_title('Shear in y direction');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("Sy [kN]");

fig4, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, Sz , linewidth=4);
ax.set_title('Shear in z direction');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("Sz [kN]");

fig5, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, T, linewidth=4);
ax.set_title('Torque');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("T [kNm]");

fig6, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, v, linewidth=4);
ax.set_title('Deflection in y direction');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("v [mm]");

fig7, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, w, linewidth=4);
ax.set_title('Deflection in z direction');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("w [mm]");

fig8, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, theta, linewidth=4);
ax.set_title('Twist');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("theta [deg]");

plt.show()
    
'''###################################################################################'''


