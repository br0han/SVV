import GlobalConstants as GC
import numpy as np
from matplotlib import pyplot as plt
from cubicspline import cubicspline


<<<<<<< HEAD
=======
data = np.genfromtxt('aerodynamicloadcrj700.dat', dtype=None, delimiter=',')
val = data[:,1]
ShearCenter = 0.2

n=len(val)
Ca = GC.Ca #0.484
la = GC.la #1.691
Coord = []

for i in range(1,n+1):
    theta0 = (i-1)/n * np.pi
    theta1 = (i)/n * np.pi
    z = -1/2*(Ca/2*(1-np.cos(theta0))+Ca/2*(1-np.cos(theta1)))
    Coord.append(z)
    
Coord = np.array(Coord)

>>>>>>> 9c13501ddf095838e1e6df9f27ff586170a31694
def DistForceTorqueMatrix():
    ###Some inputs
    data = np.genfromtxt('aerodynamicloadcrj700.dat', dtype=None, delimiter=',')
    val = data[:,1]
    ShearCenter = 0.2

    n=len(val)
    Ca = GC.Ca #0.484
    la = GC.la #1.691
    Coord = []
    for i in range(1,n+1):
        theta0 = (i-1)/n * np.pi
        theta1 = (i)/n * np.pi
        z = -1/2*(Ca/2*(1-np.cos(theta0))+Ca/2*(1-np.cos(theta1)))
        Coord.append(z)
    Coord = np.array(Coord)
    ###End of inputs


    DistributedForce = []
    DistributedTorque = []

    for i in range(0,41):   #41 because its the width of the 'data' matrix
        #----------------Part to find the point force
        val = data[:,i]
        CoefM = cubicspline(val, Coord)
        
        CoefInt = np.zeros(np.shape(CoefM))
        CoefInt[:,0] = CoefM[:,0]/4
        CoefInt[:,1] = CoefM[:,1]/3
        CoefInt[:,2] = CoefM[:,2]/2
        CoefInt[:,3] = CoefM[:,3]/1

        Force = 0
        for j in range(len(val)-1):
            Force =(CoefInt[j,0]*abs(Coord[j+1])**4 + CoefInt[j,1]*abs(Coord[j+1])**3 +CoefInt[j,2]*abs(Coord[j+1])**2+CoefInt[j,3]*abs(Coord[j+1])**1)-(CoefInt[j,0]*abs(Coord[j])**4 + CoefInt[j,1]*abs(Coord[j])**3 +CoefInt[j,2]*abs(Coord[j])**2+CoefInt[j,3]*abs(Coord[j])**1) + Force
        #print(Force)
        DistributedForce.append(Force)

        #----------------Part to find the point force line of action
        Moment = 0
        for j in range(len(val)-1):
            Moment = ((CoefInt[j,0]*abs(Coord[j+1])**4 + CoefInt[j,1]*abs(Coord[j+1])**3 +CoefInt[j,2]*abs(Coord[j+1])**2+CoefInt[j,3]*abs(Coord[j+1])**1)-(CoefInt[j,0]*abs(Coord[j])**4 + CoefInt[j,1]*abs(Coord[j])**3 +CoefInt[j,2]*abs(Coord[j])**2+CoefInt[j,3]*abs(Coord[j])**1)) * (abs(Coord[j+1])-abs(Coord[j]))/2 + Moment

        POA = Moment / Force * -1
        #print(POA)
            
        #----------------Torque genertated by the point force at the shear center

        Torque = Force * (POA - ShearCenter)  ##   Shear Center is variable
        DistributedTorque.append(Torque)


    DF  = np.array(DistributedForce)    #list of Forces on each chord (41 points)
    DT = np.array(DistributedTorque)    #list of Torques on each chord (41 points)

    n = len(DT)
    
    Ca = GC.Ca #0.484
    la = GC.la #1.691
    Span = []
    for i in range(1,n+1):
        theta0 = (i-1)/n * np.pi
        theta1 = (i)/n * np.pi
        z = 1/2*(la/2*(1-np.cos(theta0))+la/2*(1-np.cos(theta1)))
        Span.append(z)
    Span = np.array(Span)
    
<<<<<<< HEAD

    MatrixF = cubicspline(DF , Span)    #interpolation of Force along the span (40 equations)
    MatrixT = cubicspline(DT , Span)    #interpolation of Torque along the span (40 equations)
=======
    MatrixF = cubicspline(DF , Span)
    MatrixT = cubicspline(DT , Span)
>>>>>>> 9c13501ddf095838e1e6df9f27ff586170a31694


    BotRow = np.array([0,0,0,0])
    MatrixMF = np.vstack((MatrixF,BotRow))
    MatrixMT = np.vstack((MatrixT,BotRow))

    Span  = np.reshape(Span,(len(Span),1))

    MatrixMF = np.hstack((MatrixMF,Span))
    MatrixMT = np.hstack((MatrixMT,Span))
    



    

    return MatrixF, MatrixT , MatrixMF ,MatrixMT 

MatrixF, MatrixT ,MatrixMF ,MatrixMT = DistForceTorqueMatrix()
    


