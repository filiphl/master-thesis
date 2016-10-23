import re
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from os import listdir
from sys import argv



class Forces:
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
                    #self.matrix         = np.zeros([self.snChunks, self.snChunks, 3])

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
                #plt.plot(21.5,22.5, 'g*', markersize=10)
                plt.draw()

#------------------------------------------------------------------------------#




if __name__=='__main__':

    singleObject = Forces('../forceFiles/forcesAll.txt')
    singleObject.name = 'whatever'
    singleObject.plotAverage = True
    singleObject.plotMatrix()
    plt.show()
