# plot a histogram of number of particles

from matplotlib import rc
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import time
import sys
import re



if len(sys.argv) > 1:
    filepath = str(sys.argv[1])

infile   = open( filepath, "r")
content  = infile.read()
print "Measuring file size...\n"
nSamples = len(re.findall("\s1\s", content))
infile.close()
infile   = open( filepath, "r")

for i in xrange(3):
    infile.readline()

timestep, nChunks, nParticles = infile.readline().split()
snChunks = int(np.sqrt(float(nChunks))) -1
matrix = np.zeros([snChunks, snChunks])



toolbar_width = 100
<<<<<<< HEAD
ts = nSamples/100 + 1
=======
ts = nSamples/100.
>>>>>>> a7eaec87a21c68b04bad6a73a811f5ea0adc0756

sys.stdout.write("Reading data...")


# setup toolbar
sys.stdout.write("|%s|" % (" " * toolbar_width))
sys.stdout.flush()
sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['





c = 1
binWidth = 0
first = True
for line in infile:
    if first==True:
        try:
            col = line.split()
            if binWidth == 0:
                binWidth = float(col[1])*2
            x = int(float(col[1])/binWidth)
            y = int(float(col[2])/binWidth)
            matrix[x,y] = float(col[4])
        except:
            first = False
            if not c%ts:
                sys.stdout.write(u"\u25A1")
                sys.stdout.flush()
            continue
    else:
        try:
            col = line.split()
            x = int(float(col[1])/binWidth)
            y = int(float(col[2])/binWidth)
            matrix[x,y] += float(col[4])
        except:
            c+=1
            if not c%ts:
                sys.stdout.write(u"\u2588")
                sys.stdout.flush()
            continue

sys.stdout.write("\n")

print "Number of samples: ", c
print "Bin width: ", binWidth
print "Size of matrix: %dx%d"%(snChunks, snChunks)



infile.close()
matrix/=c
#matrix /= sum(sum(matrix))

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

fig = plt.figure()
plt.imshow(matrix, interpolation='none', cmap="hot_r", )
cb = plt.colorbar(format=ticker.FuncFormatter(fmt))
cb.ax.invert_yaxis()



#---------------------------- radial distribution of forces ------------------
N = 31
bins   = np.linspace(1,N, N)
values = np.zeros(N)
counts = np.zeros(N)
xc = snChunks/2.
yc = snChunks/2.

for i in xrange(np.shape(matrix)[0]):
    for j in xrange(np.shape(matrix)[1]):
        r = np.sqrt((i-xc)**2+(j-yc)**2)
        for k in xrange(len(bins)):
            if r <= bins[k]:
                values[k] += matrix[i,j]
                counts[k] += 1
                break


for i in xrange(len(counts)):
    if counts[i] != 0:
        values[i]/=counts[i]

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure()
plt.plot(bins,values, "-*", color="#8080ff", linewidth=2)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.grid("on")
plt.xlabel(r"Radial distance [\AA]")
plt.ylabel(r"Stress in z-direction [eV/\AA]")

plt.show()
