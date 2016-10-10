import re
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from os import listdir

class forces:

    def __init__(self, path, name):
        self.name = name

        infile = open(path, 'r')
        for i in xrange(3):
            infile.readline()
        self.content = infile.read()

        self.content = self.content.split('\n')
        for line in xrange(len(self.content)):
            self.content[line] = re.sub(r'\s*\d+\s\d*\.?\d*\s\d*\.?\d*\s\d+\s0\s0\s0\s?', r'', self.content[line])

        self.content = filter(None, self.content)

        self.nChunks = int(self.content[0].split()[1])
        self.snChunks = int(np.sqrt(self.nChunks))

        self.absoluteForces = np.zeros([self.snChunks, self.snChunks])
        self.matrix         = np.zeros([self.snChunks, self.snChunks, 3])
        #self.count          = np.zeros([self.snChunks, self.snChunks])

        self.binWidth = 7.121

        self.loadMatrix()



    def loadMatrix(self):
        steps = 0
        for line in self.content:
            col = line.split()
            if len(col)>3:
                x = int(round(float(col[1])/self.binWidth))
                y = int(round(float(col[2])/self.binWidth))
                fx = float(col[4])
                fy = float(col[5])
                fz = float(col[6])
                self.absoluteForces[x,y] += np.sqrt(fx**2 + fy**2 + fz**2)
            else:
                steps += 1

        self.matrix /= steps
        self.absoluteForces /= steps

class Plotter:

    def __init__(self, objects):
        if type(objects) == list:
            self.name = 'Averaged'
            self.snChunks = objects[0].snChunks
            self.absoluteForces = np.zeros([self.snChunks, self.snChunks])
            self.plotAverage(objects)

        else:
            self.name = objects.name
            self.snChunks = objects.snChunks
            self.absoluteForces = objects.absoluteForces
            self.plotMatrix()

    def plotAverage(self, objects):
        for forceObject in objects:
            self.absoluteForces += forceObject.absoluteForces
        self.absoluteForces /= len(objects)
        self.plotMatrix()

    def plotEach(self, objects):
        for forceObject in objects:
            self.absoluteForces = forceObject.absoluteForces
            self.name = forceObject.name
            self.plotMatrix()

    def fmt(self, x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    def plotMatrix(self):
        plt.imshow(self.absoluteForces, interpolation='nearest', cmap="hot_r", )
        cb = plt.colorbar(format=ticker.FuncFormatter(self.fmt))
        #cb.ax.invert_yaxis()
        plt.title(self.name)
        plt.show()

singleObject = forces('forceFiles/forces1.txt', 'Single frame t=0')
Plotter(singleObject)

singleObject = forces('forceFiles/forces9.txt', 'Single frame t=900')
Plotter(singleObject)

myObjects = []
for files in listdir('forceFiles/'):
    myObjects.append(forces('forceFiles/'+files, 'Averaged'))
Plotter(myObjects)
