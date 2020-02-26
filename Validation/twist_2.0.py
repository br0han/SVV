import pandas as pd
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt 

data = pd.read_csv('Leading_edge.txt', header = None)
data.columns = ["x","y","z"]

Leading_edge = data

X = np.array(Leading_edge["x"])

data = pd.read_csv('Tralling_edge.txt', header = None)
data.columns = ["x","y","z"]

Tralling_edge = data

u = np.array([0,1])

v = []

#calculate the vector

l = 0

while l < 108:
    
    y_1 = Leading_edge.loc[l,"y"]
    y_2 = Tralling_edge.loc[l,"y"]
    z_1 = Leading_edge.loc[l,"z"]
    z_2 = Tralling_edge.loc[l,"z"]
    dy_dz = (y_1-y_2)/(z_1-z_2)
    v = np.append(v,dy_dz)
    v = np.append(v,1)
    
    l = l+1
    
def split(array,step):
    
    for i in range(0,len(array),step):
        yield array[i:i + step]
        
def mag(x): 
    return sqrt(sum(i**2 for i in x))

v = list(split(v,2))

np.asarray(v)

#Calculate the angle of twist

mag_u = mag(u)

fun = []

j = 0

while j < 108:
    
    dot = np.dot(u,v[j])
    mag_v = mag(v[j])
    value = dot/(mag_u*mag_v)
    fun = np.append(fun,value)
    
    j = j+1
    
final_r = np.arccos(fun)
final_d = np.rad2deg(final_r)

min_deg = np.min(final_d)

twist = np.subtract(final_d,min_deg)

twist_list = np.array((X,twist))
twist_list = twist_list.transpose()
twist_list = twist_list[twist_list[:, 0].argsort()]
  
plt.plot(twist_list[:,0],twist_list[:,1]) 

plt.xlabel('X location') 
 
plt.ylabel('Twist [deg]')  
 
plt.show() 