import re
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from os import listdir
from sys import argv



class forces:

    def __init__(self, path):

        self.openFile(path)
        self.setup()
        self.setFlags()
        self.loadMatrix()


    def openFile(self, path):
        infile = open(path, 'r')
        for i in xrange(3):
            infile.readline()
        self.content = infile.read()
        self.content = self.content.split('\n')
        for line in xrange(len(self.content)):
            self.content[line] = re.sub(r'\s*\d+\s\d*\.?\d*\s\d*\.?\d*\s\d+\s0\s0\s0\s?', r'', self.content[line])
        self.content = filter(None, self.content) #Remove blank gaps


    def setup(self):
        self.nChunks        = int(self.content[0].split()[1])
        self.snChunks       = int(np.sqrt(self.nChunks))
        self.absoluteForces = np.zeros([self.snChunks, self.snChunks])
        self.matrix         = np.zeros([self.snChunks, self.snChunks, 3])
        self.binWidth       = 7.121
        self.nCells         = 0
        self.steps          = 0


    def setFlags(self):
        a = ' '.join(argv)

        self.plotAverage  = ("average".upper() in a.upper()) or ("avg".upper() in a.upper())
        self.plotEach     = ("every".upper() in a.upper()) or ("each".upper() in a.upper())
        self.save         =  "save".upper() in a.upper()
        self.animate      =  "ani".upper() in a.upper()
        self.deleteFrames =  False


        if self.animate:
            if not self.save:
                self.save = True
                self.deleteFrames = True
            if not self.plotEach:
                self.animate = False
                print "Cannot make animation without the 'each' argument."


    def loadMatrix(self):

        if self.plotEach:
            plt.figure()
            plt.ion()

        for line in self.content:
            col = line.split()

            if len(col)>3:
                x = int(round(float(col[1])/self.binWidth))
                y = int(round(float(col[2])/self.binWidth))
                fx = float(col[4])
                fy = float(col[5])
                fz = float(col[6])
                self.absoluteForces[x,y] += np.sqrt(fx**2 + fy**2 + fz**2)
                self.matrix[x,y] += [fx, fy, fz]

            else:
                if self.plotEach:
                    self.name='time step %03d'%self.steps
                    self.plotMatrix()
                    self.absoluteForces = np.zeros([self.snChunks, self.snChunks])
                    self.matrix         = np.zeros([self.snChunks, self.snChunks, 3])

                self.steps += 1
                self.nCells = int(col[1])

        self.absoluteForces /= self.steps
        self.matrix         /= self.steps

        if self.plotAverage:
            self.name='Average'
            self.plotMatrix()

        if self.animate:
            from scitools.std import movie
            print '\n\nProducing gif from image files...'
            movie('timestep*.png', encoder='convert', fps=25, output_file='forceDistribution.gif')	# Make a gif
            if self.deleteFrames:
                import os
                os.system('rm timestep*.png')
                print 'Deleted image files.'


    def fmt(self, x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)


    def plotMatrix(self):
        plt.imshow(self.absoluteForces, interpolation='nearest', cmap="hot_r", )
        plt.colorbar(format=ticker.FuncFormatter(self.fmt))
        #plt.clim(0,0.05)
        #cb = plt.colorbar(format=ticker.FuncFormatter(self.fmt))
        #cb.ax.invert_yaxis()
        plt.title(self.name)
        if self.save:
            name = "%s.png"%"".join(self.name.split())
            plt.savefig(name, format="png")
            plt.close()
            print "Saved as:", name

        else:
            if self.plotEach:
                plt.draw()
                plt.clf()

            if self.plotAverage:
                plt.show()

#------------------------------------------------------------------------------#
class surface:

    def __init__(self, path):

        self.N = 4050
        self.loadAll(path)


    def loadAll(path):
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


            p = np.zeros([self.N, 3]) #nAtoms*3
            i=0
            for line in infile:
                col = line.split()
                p[i,0] = float(col[2])
                p[i,1] = float(col[3])
                p[i,2] = float(col[4])
                i+=1

            x = np.linspace(0, 1, N+1)
            y = np.linspace(0, 1, N+1)

            for i in xrange(np.shape(p)[0]):
                for j in xrange(len(x)-1, -1, -1):
                    if p[i,0] > x[j]:
                        for k in xrange(len(y)-1, -1, -1):
                            if p[i,1] > y[k]:
                                grid[j][k].append(p[i])
                                break
                        break
        return grid






    def leastSquarePlane(points):
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

        #normalVector = Norm([1, parameter[0], parameter[1]])   # ai, bj, xk
        return parameter


    def getAngle(v1, v2):
        'Computes angle between a vector and a line parallel to another vector'
        lv1 = np.linalg.norm(v1)
        lv2 = np.linalg.norm(v2)

        angle = np.arccos(np.dot(v1,v2)/(lv1*lv2))
        if np.absolute(angle-np.pi/2)<angle:
            angle = np.absolute(angle-np.pi/2)
        return angle #[rad]


    def Norm(vector):
        for i in xrange(3):
            vector[i]/=np.linalg.norm(vector)
        return vector



    def plotCell(i,j, parameter):
        "Plot this shit"
        X = [0, 1]
        Y = [0, 1]
        X,Y = np.meshgrid(X,Y)
        Z = X*parameter[0] + Y*parameter[1] + parameter[2]
        X = [i, i+1]
        Y = [j, j+1]
        X,Y = np.meshgrid(X,Y)
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap="gray", vmin=0, vmax=1, antialiased=False)
        plt.hold("on")


    #------------------------------------------------------------------------------#



    if __name__=='__main__':

        singleObject = forces('forceFiles/forcesAll.txt')
