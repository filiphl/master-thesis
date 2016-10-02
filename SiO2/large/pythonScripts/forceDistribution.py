import re
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import sys

filepath = str(sys.argv[1])
infile = open(filepath, 'r')

for i in xrange(3):
    infile.readline()
content = infile.read()
content = content.split('\n')
for line in xrange(len(content)):
    content[line] = re.sub(r'\s*\d+\s\d*\.?\d*\s\d*\.?\d*\s\d+\s0\s0\s0\s?', r'', content[line])


content = filter(None, content)

nChunks = int(content[0].split()[1])
snChunks = int(np.sqrt(nChunks))


NormalForces = np.zeros([snChunks, snChunks])
Force = np.zeros([snChunks, snChunks, 3])
count  = np.zeros([snChunks,snChunks])

binWidth = 7.12
steps = 1
for line in content:
    col = line.split()
    if len(col)>3:
        x = int(round(float(col[1])/binWidth))
        y = int(round(float(col[2])/binWidth))
        fx = float(col[4])
        fy = float(col[5])
        fz = float(col[6])
        Force[x,y] += [fx,fy,fz]
        count [x,y] += 1
    else:
        steps += 1

ts = snChunks/100

for x in xrange(snChunks):
    for y in xrange(snChunks):
        if count[x,y] > 0:
            Force[x,y] /= count[x,y]*steps
            NormalForces[x,y] = np.linalg.norm(Force[x,y]) * np.cos( getAngle(Force[x,y], normalVector[x,y]) )
            #absoluteForces[x,y] = sum(Force[x,y])
        else:
            Force[x,y] = [0,0,0]

#absoluteForces = np.sqrt(absoluteForces)

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

plt.imshow(absoluteForces, interpolation='nearest', cmap="hot_r", )
cb = plt.colorbar(format=ticker.FuncFormatter(fmt))
#cb.ax.invert_yaxis()
plt.show()
