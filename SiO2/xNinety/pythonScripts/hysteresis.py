from forceDistribution import *
from os import listdir

folderPath = '../forceFiles/partitions/'

directory = listdir(folderPath)


#print directory

for f in directory:
    t = int(f[6:-4])
    print t
    r = 7
    N = 46
    surfN = 35
    cx=22.5
    cy=22.5
    nn=8
    bw=1
    dist = ForceDistribution(N, surfN, nn, bw, cx, cy, timeStep=t)
    dist.computeDistributions()
    dist.plotDistributions()
    plt.savefig('timeStep%06d.pdf'%t)
