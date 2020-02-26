import pandas as pd
import numpy as np

data = pd.read_csv('input_nodes.txt', header = None)
data.columns = ["node","x_n","y_n","z_n"]

nodes = data

max_value = nodes['z_n'].max()
min_value = nodes['z_n'].min()

leading_edge = nodes.loc[nodes['z_n'] == max_value]
tralling_edge = nodes.loc[nodes['z_n'] == min_value]

leading_edge_nodes = leading_edge["node"]
leading_edge_nodes = leading_edge_nodes.sub(1)
leading_edge_nodes = np.array(leading_edge_nodes)
tralling_edge_nodes = tralling_edge["node"]
tralling_edge_nodes = tralling_edge_nodes.sub(1)
tralling_edge_nodes = np.array(tralling_edge_nodes)

data = pd.read_csv('displacement_jambent.txt', delim_whitespace = True, header = None)
data.columns = ["node","magnitude","dx","dy","dz"]

change = np.array(data)

dx = change[:,2]
dy = change[:,3]
dz = change[:,4]

i = 0

X = []
Y = []
Z = []

X_1 = nodes["x_n"]
Y_1 = nodes["y_n"]
Z_1 = nodes["z_n"]

while i < 6588:
    
    x = X_1[i]+dx[i]
    X = np.append(X,x)
    y = Y_1[i]+dy[i]
    Y = np.append(Y,y)
    z = Z_1[i]+dz[i]
    Z = np.append(Z,z)
    
    i = i+1

leading_edge_new_nodes = []

j = 0

while j < 108:
    
    leading_edge_x_value = X[leading_edge_nodes[j]] 
    leading_edge_new_nodes = np.append(leading_edge_new_nodes,leading_edge_x_value)
    leading_edge_y_value = Y[leading_edge_nodes[j]] 
    leading_edge_new_nodes = np.append(leading_edge_new_nodes,leading_edge_y_value)
    leading_edge_z_value = Z[leading_edge_nodes[j]] 
    leading_edge_new_nodes = np.append(leading_edge_new_nodes,leading_edge_z_value)
    
    j = j+1
    
def split(array,step):
    
    for i in range(0,len(array),step):
        yield array[i:i + step]

Leading_edge = list(split(leading_edge_new_nodes,3))

Leading_edge = np.asarray(Leading_edge)

tralling_edge_new_nodes = []

k = 0

while k < 108:
    
    tralling_edge_x_value = X[tralling_edge_nodes[k]] 
    tralling_edge_new_nodes = np.append(tralling_edge_new_nodes,tralling_edge_x_value)
    tralling_edge_y_value = Y[tralling_edge_nodes[k]] 
    tralling_edge_new_nodes = np.append(tralling_edge_new_nodes,tralling_edge_y_value)
    tralling_edge_z_value = Z[tralling_edge_nodes[k]] 
    tralling_edge_new_nodes = np.append(tralling_edge_new_nodes,tralling_edge_z_value)
    
    k = k+1
    
Tralling_edge = list(split(tralling_edge_new_nodes,3))

Tralling_edge = np.asarray(Tralling_edge)

#np.savetxt("Leading_edge.txt",LE,delimiter = ",")
#np.savetxt("Tralling_edge.txt",TE,delimiter = ",")

df = pd.DataFrame (Leading_edge)

filepath = 'Leading_edge.xlsx'

df.to_excel(filepath, index=False)

df = pd.DataFrame (Tralling_edge)

filepath = 'Tralling_edge.xlsx'

df.to_excel(filepath, index=False)



