from surface import *
from compareMatrices import *
from radialBinning2 import smooth
from matplotlib import rc
import pickle



def loadSurface(filePath=False, N=45, s=5):
    if filePath:
        filePath = filePath.rstrip('.pkl') + 'N%ds%d.pkl'%(N,s)
        print filePath
        try:
            with open(filePath, 'rb') as input:
                print "Loaded surface file with N=%d and s=%d."%(N,s)
                return pickle.load(input)
        except:
            print "Couldn't load surface file."
    s = SurfaceRegression('../dumpFiles/', N, False, s)
    if filePath:
        print filePath
        with open(filePath, 'wb') as output:
            pickle.dump(s, output, pickle.HIGHEST_PROTOCOL)
    return s
#surf.plotPlanes()

def loadForces(filePath=False):
    if filePath:
        try:
            with open(filePath, 'rb') as input:
                print "loaded force file."
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


#------------------------------------------------------------------------------#
N = 45
surfN = 35
nn = 5
mapping = float(N)/surfN

surf = loadSurface('../dataFiles/surface.pkl', N=surfN, s=nn)
force = loadForces('../dataFiles/forces.pkl')
radialBinning = smooth(N, (N-1)/2., binWidth=1.2)
print np.shape(force.absoluteForces), np.shape(radialBinning.weights)
#surf.plotPlanes()
#plt.figure()
#force.plotMatrix()
#plt.show()
#radialBinning.show(3)


for i in xrange(45):
    for j in xrange(45):
        #if force.matrix[i][j].any():
        #print int(i/mapping), int(j/mapping), np.cos(surf.getAngle( force.matrix[i][j], surf.grid[int(i/mapping)][int(j/mapping)] )) , "      ", force.matrix[i][j], "      " , surf.grid[int(i/mapping)][int(j/mapping)]
        if not np.isnan( np.cos( surf.getAngle( force.matrix[i][j], surf.grid[int(i/mapping)][int(j/mapping)] ) ) ):
            force.absoluteForces[i][j] *= np.cos( surf.getAngle( force.matrix[i][j], surf.grid[int(i/mapping)][int(j/mapping)] ) )
            print np.cos( surf.getAngle( force.matrix[i][j], surf.grid[int(i/mapping)][int(j/mapping)] ) )
#force.plotMatrix()


plt.figure()
force.plotMatrix()



r = radialBinning.bins[1:]
F = np.zeros(radialBinning.nBins)

for k in xrange(radialBinning.nBins):
    area = np.pi * (radialBinning.bins[k+1]**2-radialBinning.bins[k]**2)
    value = 0
    for i in xrange(45):
        for j in xrange(45):
            if not np.isnan(force.absoluteForces[i,j]):
                value += force.absoluteForces[i,j] * radialBinning.weights[i,j,k]
                #print force.absoluteForces[i,j], radialBinning.weights[i,j,k]
    F[k] = value/area

print F

#------------------------------------------------------------------------------#



rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure()
plt.plot(r,F,
 '-h',
 color="#478684",
 linewidth=2,
 markersize=7,
 markerfacecolor='#3b617c',
 fillstyle='full')
plt.title(r"$ $Radial distribution of normal force", fontsize=18)
plt.xlabel(r"$r$", fontsize=18)
plt.ylabel(r"$eV/\AA$", fontsize=16)
plt.xlim([0,12 ])
plt.grid('on')
plt.show()
