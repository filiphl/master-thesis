import numpy as np
import matplotlib.pyplot as plt
import sys

filepath = str(sys.argv[1])
infile = open(filepath, 'r')

for i in xrange(9):
    infile.readline()


p = np.zeros([4050, 3])
i=0
for line in infile:
    col = line.split()
    p[i,0] = float(col[2])
    p[i,1] = float(col[3])
    p[i,2] = float(col[4])
    i+=1

N = 22
grid = []
for i in xrange(N):
    grid.append([])
    for j in xrange(N):
        grid[i].append([])

x = np.linspace(0,1, N+1)
y = np.linspace(0,1, N+1)

for i in xrange(np.shape(p)[0]):
    for j in xrange(len(x)-1, -1, -1):
        if p[i,0] > x[j]:
            for k in xrange(len(y)-1, -1, -1):
                if p[i,1] > y[k]:
                    grid[j][k].append(int(i))
                    break
            break


for row in grid:
    for cell in row:
        print len(cell),

def leastSquarePlane(points):
    """ Function that uses least square linear regresion to approximate a plane
    from a set of (x,y,z) coordinates."""

    n = len(points)

    X  = 0;     Y  = 0;     Z  = 0
    X2 = 0;     XY = 0;     XZ = 0
    Y2 = 0;     YZ = 0

    for point in points:
        x = point[0]
        y = point[1]
        z = point[2]
        X += x
        Y += y
        Z += z
        XX += x*x
        XY += x*y
        XZ += x*z
        YY += y*y
        YZ += y*z

    A = np.matrix([[X2, XY, X],[XY, Y2, Y],[X, Y, n]])
    b = np.array([XZ, YZ, Z])




    return a,b,c
