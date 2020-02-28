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
    T[i] = -graph.T(x[i])
    v[i] = graph.v(x[i])
    w[i] = graph.w(x[i])
    theta[i] = -graph.theta(x[i])*(180/np.pi)
    
print ("That took", time.time() - t0, "second(s)")

# plt.plot(x, theta)
# plt.grid(1)
# plt.show()
    
    
fig, ax = plt.subplots(1, 3, figsize = (24, 8));
ax[0].plot(x, My);
ax[0].set_title('Moment in y direction');
ax[0].grid(True);
ax[0].set_xlabel("x [m]");
ax[0].set_ylabel("My [Nm]");

ax[1].plot(x, w);
ax[1].set_title('Deflection in z direction');
ax[1].grid(True);
ax[1].set_xlabel("x [m]");
ax[1].set_ylabel("w [m]");

ax[2].plot(x, Sz);
ax[2].set_title('Shear force in z direction');
ax[2].grid(True);
ax[2].set_xlabel("x [m]");
ax[2].set_ylabel("Sz [N]");
plt.tight_layout()



fig, ax = plt.subplots(1, 3, figsize = (24, 8));
ax[0].plot(x, Mz);
ax[0].set_title('Moment in z direction');
ax[0].grid(True);
ax[0].set_xlabel("x [m]");
ax[0].set_ylabel("Mz [Nm]");

ax[1].plot(x, v);
ax[1].set_title('Deflection in y direction');
ax[1].grid(True);
ax[1].set_xlabel("x [m]");
ax[1].set_ylabel("v [m]");

ax[2].plot(x, Sy);
ax[2].set_title('Shear force in y direction');
ax[2].grid(True);
ax[2].set_xlabel("x [m]");
ax[2].set_ylabel("Sz [N]");
plt.tight_layout()



fig, ax = plt.subplots(1, 2, figsize = (24, 8));
ax[0].plot(x, T);
ax[0].set_title('Torque');
ax[0].grid(True);
ax[0].set_xlabel("x [m]");
ax[0].set_ylabel("T [Nm]");

ax[1].plot(x, theta);
ax[1].set_title('Twist');
ax[1].grid(True);
ax[1].set_xlabel("x [m]");
ax[1].set_ylabel("theta [m]");
plt.tight_layout()




