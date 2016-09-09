from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np




# -----------------------------------------------------------------------------

# plot a histogram of number of particles

from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from sys import argv

group = 1
if len(argv) > 1:
    group = float(argv[1])

infile = open("forcesInChunk.txt", "r")
for i in xrange(3):
    print infile.readline()
timestep, nChunks, nParticles = infile.readline().split()

snChunks = int(np.sqrt(float(nChunks)))/group

matrix = np.zeros([snChunks, snChunks])



c = 0
binWidth = 0
first = True
for line in infile:
    if first==True:
        try:
            col = line.split()
            if binWidth == 0:
                binWidth = float(col[1])*2
            x = int(float(col[1])/(binWidth*group))
            y = int(float(col[2])/(binWidth*group))
            matrix[x,y] = float(col[4])
        except:
            first = False
            c+=1
            continue
    else:
        try:
            col = line.split()
            x = int(float(col[1])/(binWidth*group))
            y = int(float(col[2])/(binWidth*group))
            matrix[x,y] += float(col[4])
        except:
            c+=1
            continue

print "Number of samples: ", c
print "Bin width: ", binWidth
infile.close()
matrix/=c
#matrix /= sum(sum(matrix))

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


plt.imshow(matrix, interpolation='none', cmap="hot_r", )
cb = plt.colorbar(format=ticker.FuncFormatter(fmt))
cb.ax.invert_yaxis()
plt.show()

#----------------------------------------------------------------------------

X = np.linspace(0, snChunks, snChunks)
Y = np.linspace(0, snChunks, snChunks)

X,Y = np.meshgrid(X,Y)
print np.size(X), np.size(Y), np.size(matrix)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, matrix, rstride=1, cstride=1, cmap="hot_r",
                       linewidth=0, antialiased=False)
#ax.set_axis_bgcolor('#cccccc')
ax.w_xaxis.set_pane_color((0.8, 0.8, 0.8, 1.0))
ax.w_yaxis.set_pane_color((0.8, 0.8, 0.8, 1.0))
ax.w_zaxis.set_pane_color((0.8, 0.8, 0.8, 1.0))
#ax.set_zlim(-1.01, 1.01)

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
