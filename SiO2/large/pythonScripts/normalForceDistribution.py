from surface import *
from compareMatrices import *
from matplotlib import rc
import pickle

N = 45

mapping = 45.0/N

def loadSurface(filePath=False):
    if filePath:
        try:
            with open(filePath, 'rb') as input:
                print "loaded surface file."
                return pickle.load(input)
        except:
            print "Couldn't load surface file."
    s = SurfaceRegression('../dumpFiles/', N, False, 8)
    with open(filePath, 'wb') as output:
        return pickle.dump(s, output, pickle.HIGHEST_PROTOCOL)
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
    with open(filePath, 'wb') as output:
        return pickle.dump(F, output, pickle.HIGHEST_PROTOCOL)
    return F


surf = loadSurface('../dataFiles/surface.pkl')

force = loadForces('../dataFiles/forces.pkl')
force.plotMatrix()

for i in xrange(45):
    for j in xrange(45):
        #if force.matrix[i][j].any():
            #print int(i/mapping), int(j/mapping), np.cos(surf.getAngle( force.matrix[i][j], surf.grid[int(i/mapping)][int(j/mapping)] )) , "      ", force.matrix[i][j], "      " , surf.grid[int(i/mapping)][int(j/mapping)]
        force.absoluteForces[i][j] *= np.cos( surf.getAngle( force.matrix[i][j], surf.grid[int(i/mapping)][int(j/mapping)] ) )

#force.plotMatrix()



xc = 22.5
yc = 22.5
R  = 10
M  = 10

F     = np.zeros(M)
count = np.zeros(M)

r = np.linspace(1, R, M)

for i in xrange(45):
    for j in xrange(45):
        cr = np.sqrt((i-xc)**2 + (j-yc)**2)
        for k in xrange(M):
            if cr < r[k]:
                count[k]+=1
                break
        if cr < R:
            F[k]+=force.absoluteForces[i][j]

for i in xrange(M):
    if count[i]:
        F[i]/=count[i]
        if not F[i] > 0:
            F[i] = 0
#------------------------------------------------------------------------------#



rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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
plt.xlim([0,R])
plt.grid('on')
plt.show()