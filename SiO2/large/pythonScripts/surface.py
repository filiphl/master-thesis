import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import sys
from os import listdir


class surfaceRegression:

    def __init__(self, path):

        self.N = 15
        self.numberOfAtoms = 4050

        self.loadSurfacePoints(path)
        self.approximatePlanes()
        self.plotPlanes()





    def loadSurfacePoints(self, path):
        self.grid = []
        for i in xrange(self.N):
            self.grid.append([])
            for j in xrange(self.N):
                self.grid[i].append([])

        directory = listdir(path)
        for files in directory:
            infile = open(path+files, 'r')

            for i in xrange(9): # dump file
                infile.readline()


            p = np.zeros([self.numberOfAtoms, 3]) #nAtoms*3
            i=0
            for line in infile:
                col = line.split()
                if float(col[4]) < 54:
                    p[i,0] = float(col[2])
                    p[i,1] = float(col[3])
                    p[i,2] = float(col[4])
                    i+=1

            x = np.linspace(0, 1, self.N+1)
            y = np.linspace(0, 1, self.N+1)

            for i in xrange(np.shape(p)[0]):
                for j in xrange(len(x)-1, -1, -1):
                    if p[i,0] > x[j]:
                        for k in xrange(len(y)-1, -1, -1):
                            if p[i,1] > y[k]:
                                self.grid[j][k].append(p[i])
                                break
                        break


    def approximatePlanes(self):
        for i in xrange(self.N):
            for j in xrange(self.N):
                #print len(self.grid[i][j])
                self.grid[i][j] = self.leastSquarePlane(self.grid[i][j])

        self.grid = np.asarray(self.grid)




    def leastSquarePlane(self, points):
        """ Function that uses least square linear regresion to approximate a plane
        from a set of (x,y,z) coordinates."""

        n = len(points)

        X  = 0;     Y  = 0;     Z  = 0
        XX = 0;     XY = 0;     XZ = 0
        YY = 0;     YZ = 0

        for point in points:
            x = point[0]
            y = point[1]
            z = point[2]
            X += x
            Y += y
            Z += z
            XX += x*x
            XY += x*y
            XZ += x*z
            YY += y*y
            YZ += y*z

        A = np.matrix([[XX, XY, X], [XY, YY, Y], [X, Y, n]])
        b = np.array([XZ, YZ, Z])

        parameter = np.linalg.solve(A,b)

        #normalVector = Norm([1, self.parameter[0], self.parameter[1]])   # ai, bj, xk
        return parameter


    def getAngle(self, v1, v2):
        'Computes angle between a vector and a line parallel to another vector'
        lv1 = np.linalg.norm(v1)
        lv2 = np.linalg.norm(v2)

        angle = np.arccos(np.dot(v1,v2)/(lv1*lv2))
        if np.absolute(angle-np.pi/2)<angle:
            angle = np.absolute(angle-np.pi/2)
        return angle #[rad]


    def Norm(self, vector):
        for i in xrange(3):
            vector[i]/=np.linalg.norm(vector)
        return vector



    def plotPlanes(self):

        fig = plt.figure()
        ax  = fig.gca(projection='3d')

        for i in xrange(self.N):
            print "i =", i
            for j in xrange(self.N):
                "Plot this shit"
                X = [0, 1]
                Y = [0, 1]
                X,Y = np.meshgrid(X,Y)
                Z = X*self.grid[i][j][0] + Y*self.grid[i][j][1] + self.grid[i][j][2]
                X = [i, i+1]
                Y = [j, j+1]
                X,Y = np.meshgrid(X,Y)
                surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap="gray", vmin=0, vmax=1, antialiased=False)
                plt.hold("on")

        plt.xlabel('x')
        plt.ylabel('y')
        ax.set_zlim(-10, 10)
        plt.show()
# ---------------------------------------------------------------------------- #

if __name__ == '__main__':
    try:
        path = str(sys.argv[1])
    except:
        print 'Please provide a valid path as an argument.'
        exit()

    surf = surfaceRegression(path)

# ---------------------------------------------------------------------------- #


#ax.set_axis_bgcolor('#cccccc')
#ax.set_zlim(-1.01, 1.01)

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

#fig.colorbar(surf, shrink=0.5, aspect=5)
