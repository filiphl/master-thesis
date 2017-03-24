from matplotlib import pyplot as plt
import numpy as np
from sys import argv

filename = str(argv[1])


infile = open(filename, 'r')

for i in xrange(3):
    infile.readline()

N = 100 #int(infile.readline().split()[1])
Nx = Ny = int(np.sqrt(N))
L = 72.1*2

stress = []
for i in xrange(Nx):
    stress.append([])
    for j in xrange(Ny):
        stress[i].append([])

x = np.linspace(0, L, Nx+1)
y = np.linspace(0, L, Nx+1)

for line in infile:
    col = line.split()
    cx = float(col[1])
    cy = float(col[2])
    for i in xrange(Nx+1):
        if x[i]>cx:
            for j in xrange(Ny+1):
                if y[j]>cy:
                    print col
                    stress[i-1][j-1].append(float(col[-1]))
                    break
            break
for i in xrange(Nx):
    for j in xrange(Ny):
        stress[i][j]=np.mean(stress[i][j])

print np.min(stress)
print np.max(stress)
stress = stress/np.min(stress)


plt.imshow(stress, interpolation='none', cmap="hot_r")
cb = plt.colorbar()
#cb.ax.invert_yaxis()
plt.show()
