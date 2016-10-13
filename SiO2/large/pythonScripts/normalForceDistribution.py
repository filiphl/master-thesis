from surface import *
from compareMatrices import *

surf = SurfaceRegression('../dumpFiles/', 45, False, 8)

force = Forces('../forceFiles/forcesAll.txt')
force.name = 'name'
force.plotAverage = True

#force.plotMatrix()

for i in xrange(45):
    for j in xrange(45):
        force.absoluteForces[i][j] *= np.cos( surf.getAngle( force.matrix[i][j], surf.grid[i][j] ) )

force.plotMatrix()
