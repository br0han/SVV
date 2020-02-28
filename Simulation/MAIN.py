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


# plt.figure(1)
# plt.plot(x, My);
# plt.set_title('Moment in y direction');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("My [Nm]");

# plt.figure(1)
# plt.plot(x, Mz);
# plt.set_title('Moment in z direction');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("Mz [Nm]");

# plt.figure(1)
# plt.plot(x, Sy);
# plt.set_title('Shear in y direction');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("Sy [Nm]");

# plt.figure(1)
# plt.plot(x, Sz);
# plt.set_title('Shear in z direction');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("Sz [Nm]");

# plt.figure(1)
# plt.plot(x, T);
# plt.set_title('Torque');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("T [Nm]");

# plt.figure(1)
# plt.plot(x, v);
# plt.set_title('Deflection in y direction');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("v [Nm]");

# plt.figure(1)
# plt.plot(x, w);
# plt.set_title('Deflection in z direction');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("w [Nm]");

# plt.figure(1)
# plt.plot(x, theta);
# plt.set_title('Twist');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("theta [Nm]");

# plt.show()
    
'''#############################################################################'''

# plt.figure(figsize = (28,0))
# plt.plot(x, My);
# plt.set_title('Moment in y direction');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("My [Nm]");


# plt.figure(figsize = (28,0))
# plt.plot(x, Mz);
# plt.set_title('Moment in z direction');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("Mz [Nm]");

# plt.figure(figsize = (28,0))
# plt.plot(x, Sy);
# plt.set_title('Shear in y direction');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("Sy [Nm]");

# plt.figure(figsize = (28,0))
# plt.plot(x, Sz);
# plt.set_title('Shear in z direction');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("Sz [Nm]");

# plt.figure(figsize = (28,0))
# plt.plot(x, T);
# plt.set_title('Torque');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("T [Nm]");

# plt.figure(figsize = (28,0))
# plt.plot(x, v);
# plt.set_title('Deflection in y direction');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("v [Nm]");

# plt.figure(figsize = (28,0))
# plt.plot(x, w);
# plt.set_title('Deflection in z direction');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("w [Nm]");

# plt.figure(figsize = (28,0))
# plt.plot(x, theta);
# plt.set_title('Twist');
# plt.grid(True);
# plt.set_xlabel("x [m]");
# plt.set_ylabel("theta [Nm]");

# plt.show()

'''#####################################################################################'''

fig1, ax1 = plt.subplots(1, figsize = (8, 8))
ax1.plot(x, My);
ax1.set_title('Moment in y direction');
ax1.grid(True);
ax1.set_xlabel("x [m]");
ax1.set_ylabel("My [kNm]");
plt.tight_layout()

fig2, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, Mz);
ax.set_title('Moment in z direction');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("Mz [kNm]");

fig3, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, Sy);
ax.set_title('Shear in y direction');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("Sy [kN]");

fig4, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, Sz);
ax.set_title('Shear in z direction');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("Sz [kN]");

fig5, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, T);
ax.set_title('Torque');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("T [kNm]");

fig6, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, v);
ax.set_title('Deflection in y direction');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("v [mm]");

fig7, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, w);
ax.set_title('Deflection in z direction');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("w [mm]");

fig8, ax = plt.subplots(1, figsize = (8, 8))
ax.plot(x, theta);
ax.set_title('Twist');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("theta [deg]");

plt.show()
    
'''###################################################################################'''
   
# fig, ax = plt.subplots(1, 3, figsize = (24, 8));
# ax[0].plot(x, My);
# ax[0].set_title('Moment in y direction');
# ax[0].grid(True);
# ax[0].set_xlabel("x [m]");
# ax[0].set_ylabel("My [Nm]");

# ax[1].plot(x, w);
# ax[1].set_title('Deflection in z direction');
# ax[1].grid(True);
# ax[1].set_xlabel("x [m]");
# ax[1].set_ylabel("w [m]");

# ax[2].plot(x, Sz);
# ax[2].set_title('Shear force in z direction');
# ax[2].grid(True);
# ax[2].set_xlabel("x [m]");
# ax[2].set_ylabel("Sz [N]");
# plt.tight_layout()



# fig, ax = plt.subplots(1, 3, figsize = (24, 8));
# ax[0].plot(x, Mz);
# ax[0].set_title('Moment in z direction');
# ax[0].grid(True);
# ax[0].set_xlabel("x [m]");
# ax[0].set_ylabel("Mz [Nm]");

# ax[1].plot(x, v);
# ax[1].set_title('Deflection in y direction');
# ax[1].grid(True);
# ax[1].set_xlabel("x [m]");
# ax[1].set_ylabel("v [m]");

# ax[2].plot(x, Sy);
# ax[2].set_title('Shear force in y direction');
# ax[2].grid(True);
# ax[2].set_xlabel("x [m]");
# ax[2].set_ylabel("Sz [N]");
# plt.tight_layout()


'''
fig, ax = plt.subplots(1, figsize = (8, 8));
ax.plot(x, T);
ax.set_title('Torque');
ax.grid(True);
ax.set_xlabel("x [m]");
ax.set_ylabel("T [Nm]");

# # ax[1].plot(x, theta);
# # ax[1].set_title('Twist');
# # ax[1].grid(True);
# # ax[1].set_xlabel("x [m]");
# # ax[1].set_ylabel("theta [m]");
plt.tight_layout()

'''


