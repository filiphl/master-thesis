from surface import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import numpy as np
import sys


def loadSurfacePoints(path, N=46):
    grid = []
    zPosition = []
    for i in xrange(N):
        grid.append([])
        zPosition.append([])
        for j in xrange(N):
            grid[i].append([])
            zPosition[i].append([])


    directory = listdir(path)
    for files in directory:
        infile = open(path+files, 'r')

        for i in xrange(9): # dump file
            infile.readline()


        p = np.zeros([4232, 3]) #nAtoms*3
        i=0
        for line in infile:
            col = line.split()
            #if float(col[4]) < 0.39:            # Only atoms at z<54
            p[i,0] = float(col[2])
            p[i,1] = float(col[3])
            p[i,2] = float(col[4])
            i+=1

        x = np.linspace(0, 1, N+1)
        y = np.linspace(0, 1, N+1)

        for i in xrange(np.shape(p)[0]):
            for j in xrange(len(x)-1, -1, -1):
                if p[i,0] > x[j]:
                    for k in xrange(len(y)-1, -1, -1):
                        if p[i,1] > y[k]:
                            grid[j][k].append(p[i])
                            zPosition[j][k].append(p[i][2])
                            break
                    break
    return grid

surface = loadSurfacePoints('../surfaceFiles/m3/Surface86000_m3/')


N = len(surface)

fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
for x in xrange(N):
    print " %d %% \r"%(x*100/N),
    sys.stdout.flush()
    m = 0
    c=0
    for y in xrange(N/2-1, N/2+2):
        for z in surface[x][y]:
            m+=z[2]
            c+=1
    if c>1:
        surface[x] = m/c
    else:
        surface[x] = m

for x in range(N):
    if surface[x] == 0:
        print "Did some shady stuff!"
        surface[x] = (surface[x-1] + surface[x+1])/2



print np.shape(surface)
plt.plot(range(N), surface)
plt.ylim([0,1])
plt.show()
