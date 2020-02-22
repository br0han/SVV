import pandas as pd
import numpy as np

data = pd.read_csv('input_nodes.txt', header = None)
data.columns = ["node","x_n","y_n","z_n"]

nodes = data

data = pd.read_csv('elements.txt', header=None)
data.columns = ["number_element", "node_1", "node_2", "node_3","node_4"]

elements = data

elements_num = elements.number_element

elements_node_1 = elements.node_1
elements_node_2 = elements.node_2
elements_node_3 = elements.node_3
elements_node_4 = elements.node_4

i = 0

element_x = []
element_y = []
element_z = []

while i < len(elements_num):

    x = []
    y = []
    z = []    
            
    n1 = elements_node_1[i]
    n2 = elements_node_2[i]
    n3 = elements_node_3[i]
    n4 = elements_node_4[i]
    
    x_loc = nodes.loc[n1-1,'x_n']
    x = np.append(x,x_loc)
    x_loc = nodes.loc[n2-1,'x_n']
    x = np.append(x,x_loc)
    x_loc = nodes.loc[n3-1,'x_n']
    x = np.append(x,x_loc)
    x_loc = nodes.loc[n4-1,'x_n']
    x = np.append(x,x_loc)
    
    y_loc = nodes.loc[n1-1,'y_n']
    y = np.append(y,y_loc)
    y_loc = nodes.loc[n2-1,'y_n']
    y = np.append(y,y_loc)
    y_loc = nodes.loc[n3-1,'y_n']
    y = np.append(y,y_loc)
    y_loc = nodes.loc[n4-1,'y_n']
    y = np.append(y,y_loc)
    
    z_loc = nodes.loc[n1-1,'z_n']
    z = np.append(z,z_loc)
    z_loc = nodes.loc[n2-1,'z_n']
    z = np.append(z,z_loc)
    z_loc = nodes.loc[n3-1,'z_n']
    z = np.append(z,z_loc)
    z_loc = nodes.loc[n4-1,'z_n']
    z = np.append(z,z_loc)       
    
    value_x = np.mean(x)
    value_y = np.mean(y)
    value_z = np.mean(z)
    
    element_x = np.append(element_x,value_x)
    element_y = np.append(element_y,value_y)
    element_z = np.append(element_z,value_z)
     
    i = i+1

i = 0

elements_list = []

while i < len(elements_num):
    
    elements_list = np.append(elements_list,elements_num[i])
    elements_list = np.append(elements_list,element_x[i])
    elements_list = np.append(elements_list,element_y[i])
    elements_list = np.append(elements_list,element_z[i])
    
    i = i+1
   
def split(array,step):
    
    for i in range(0,len(array),step):
        yield array[i:i + step]

element = list(split(elements_list,4))

np.asarray(element)

#np.genfromtxt("elements_locations",dtype=<class'float',comments='#',delimiter',',skip_header=0,skip_footer=0,c)
              
np.savetxt("elements_location.txt",element,delimiter = ",")