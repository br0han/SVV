from matplotlib import pyplot as plt
import numpy as np
# Every parameter should be unit of [m].
import GlobalConstants as G

def Coord_out(Ca = G.Ca, h = G.ha, Nst = G.nst):


    r = h / 2 # [m]

    Coord = np.zeros((Nst, 2))

    # --------- Adding the first coordinate at LE

    Coord[0, 0] = r
    Coord[0, 1] = 0

    # --------- Parameter and length between the booms

    Par = np.pi * r + 2 * np.sqrt((Ca - r) ** 2 + r ** 2)
    Dis = Par / Nst

    # --------- Check how many more is in the half circle

    Theta = Dis * 360 / (2 * np.pi * r)

    IHC = int((90 - np.mod(90, Theta)) / Theta)

    # --------- Coordinates in Circle

    for i in range(1, IHC + 1):
        Coord[i, 0] = r * np.cos(i * Theta * (2 * np.pi) / 360)
        Coord[i, 1] = r * np.sin(i * Theta * (2 * np.pi) / 360)

    for j in range(1, int(((Nst - 1) / 2 - IHC)) + 1):
        a = r / (Ca - r)
        b = r

        Diag = Dis / 2 + Dis * (j - 1)

        c = 1 + a ** 2

        lenz = np.sqrt(Diag ** 2 / c)

        z = -b / a + lenz
        y = z * a + b

        Coord[int(((Nst - 1) / 2 - j + 1)), 0] = z
        Coord[int(((Nst - 1) / 2 - j + 1)), 1] = y

    for k in range(0, int((Nst - 1) / 2)):
        Coord[Nst - 1 - k, 0] = Coord[k + 1, 0]
        Coord[Nst - 1 - k, 1] = -Coord[k + 1, 1]

    #Coord[:, 0] = Coord[:,0]- 0.0865 # Shifts to origin at leading edge, if needed.
    #Coord[:, 0] = -Coord[:, 0]
    return Coord
Coord = Coord_out()

X = Coord[:, 0]
Y = Coord[:, 1]

#plt.plot(X,Y,'o')
#plt.axis('equal')
#plt.show()

#print(sum(Coord_out(0.484)[1:3,1]))
