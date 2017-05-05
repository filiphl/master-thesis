from compareMatrices import Forces
import numpy as np
from scipy import ndimage
from matplotlib import pyplot as plt
from matplotlib import rc
import matplotlib as mpl
import pickle

def loadForces(filePath=False):
    if filePath:
        try:
            with open(filePath, 'rb') as input:
                print "Loaded force file."
                return pickle.load(input)
        except:
            print "Couldn't load force file."
    F = Forces('../forceFiles/forcesAll.txt')
    F.plotAverage = True
    F.name = 'Averaged normal force'
    if filePath:
        with open(filePath, 'wb') as output:
            pickle.dump(F, output, pickle.HIGHEST_PROTOCOL)
    return F


singleObject = loadForces(filePath='../dataFiles/forces.pkl')
singleObject.name = 'each'
myMatrix = singleObject.absoluteForces
print np.shape(myMatrix)

N = 45
M = 45
R = 15
cx = 22
cy = 22
r = np.linspace(0,R,N)
theta = np.linspace(0, 2*np.pi, M)
myevalmatrix = np.zeros((N, M, 2))
for i in range(N):
    for j in range(M):
        myevalmatrix[i, j,:] = np.asarray([cx+r[i]*np.cos(theta[j]), cy+r[i]*np.sin(theta[j])])
#output = np.zeros((M, N))
#print myevalmatrix[2, :]

output = ndimage.map_coordinates(myMatrix, np.transpose(myevalmatrix[:, :]), order=1)




#---------------------------------- Plot --------------------------------------#
fig, ax = plt.subplots(2,1, sharex=True)

im=ax[0].pcolor(output)
ax[0].axis([0,N,0,M])
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.89, 0.535, 0.02, 0.3648])
cb = fig.colorbar(im, cax=cbar_ax)

radialDist = np.mean(output, 0) #Radial mean
ax[1].plot(radialDist, linewidth=2, color="#478684")


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax[0].set_ylabel('theta', fontsize=14)
ax[1].set_ylabel('eV/A',  fontsize=14)
ax[1].set_xlabel('r',     fontsize=14)
n=6
m=5
ax[0].set_yticks(np.linspace(0,M,m))
ax[0].set_yticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'],  fontsize=14)
ax[1].set_xticks(np.linspace(0,N,n))
ax[1].set_xticklabels(['%.0f'%i for i in np.linspace(0,R,n)])
ax[1].grid('on')

cb.update_ticks()
#singleObject.plotAverage = True
#singleObject.plotMatrix()
plt.show()
