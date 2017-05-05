from matplotlib import pyplot as plt
import numpy as np

infile = open('force.txt', 'r')

for line in xrange(3):
    infile.readline()

N = int(infile.readline().split()[1])
Nx = Ny = int(np.sqrt(N))
L = 71.2

stress = []
for i in xrange(Nx):
    stress.append([])
    for j in xrange(Ny):
        stress[i].append([])
x = np.linspace(0, L, Nx+1)
y = np.linspace(0, L, Nx+1)
