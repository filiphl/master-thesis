
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import numpy
from numpy.random import randn
from scipy import array, newaxis

filepath = sys.argv[1]
infile = open(filepath, "r")

for i in xrange(10):
    infile.readline()

x = []
y = []
z = []

for line in infile:
    col = line.split()
    print col
    x.append(float(col[2]))
    y.append(float(col[3]))
    z.append(float(col[4]))
    if (z[-1]>0.44 or x[-1]>0.5):
        del x[-1]
        del y[-1]
        del z[-1]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0)
fig.colorbar(surf)

ax.xaxis.set_major_locator(MaxNLocator(5))
ax.yaxis.set_major_locator(MaxNLocator(6))
ax.zaxis.set_major_locator(MaxNLocator(5))
ax.axis('equal')
ax.set_zlim([0,1])
fig.tight_layout()
ax.view_init(20, 50)

plt.show() # or:
