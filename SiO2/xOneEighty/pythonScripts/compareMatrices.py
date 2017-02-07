import re
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from os import listdir
from sys import argv
from scipy import ndimage
import matplotlib as mpl
from matplotlib import rc

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
        self.binWidth       = 7.120
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
        for line in self.content:
            col = line.split()

            if len(col)>3:
                x = int(float(col[1])/self.binWidth)
                y = int(float(col[2])/self.binWidth)
                fx = float(col[4])
                fy = float(col[5])
                fz = float(col[6])
                self.absoluteForces[x,y] += np.sqrt(fx**2 + fy**2 + fz**2)
                self.matrix[x,y] += [fx, fy, fz]

            else:
                if self.plotEach:
                    self.name='time step %d'%int(col[0])
                    self.plotMatrix()
                    self.absoluteForces = np.zeros([self.snChunks, self.snChunks])
                    #self.matrix         = np.zeros([self.snChunks, self.snChunks, 3])

                self.steps += 1
                self.nCells = int(col[1])
        for i in xrange(np.shape(self.matrix)[0]):
            for j in xrange(np.shape(self.matrix)[1]):
                self.absoluteForces[i,j] = np.sqrt(sum(self.matrix[i,j]**2))
        self.absoluteForces /= self.steps


        self.matrix /= self.steps

        if self.plotAverage:
            self.name='Average'

        #self.plotMatrix()

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
        self.fig, self.ax = plt.subplots(3,1, figsize=(3.75,10))

        N = 100
        M = 100
        R = 16
        r = np.linspace(0,R,N)
        theta = np.linspace(0, 2*np.pi, M)
        myevalmatrix = np.zeros((N, M, 2))
        for i in range(N):
            for j in range(M):
                myevalmatrix[i, j,:] = np.asarray([22.5+r[i]*np.cos(theta[j]), 22.5+r[i]*np.sin(theta[j])])
        output = ndimage.map_coordinates(self.absoluteForces, np.transpose(myevalmatrix[:, :]), order=1)
        radialDist = np.mean(output, 0)

        im0 = self.ax[0].imshow(self.absoluteForces, interpolation='nearest', cmap="hot_r")
        im1 = self.ax[1].pcolor(output)
        self.ax[2].plot(radialDist, linewidth=2, color="#478684")

        a = 100
        p0 = max(radialDist)
        pd = lambda p0, r, a: p0*np.sqrt(1-np.linspace(0,a,100)**2/a**2)
        p  = pd(p0, r, a)

        self.ax[2].hold('on')
        self.ax[2].plot(np.linspace(0,100,100), p)
        #self.ax[1].axis([0,N,0,M])

        #self.fig.subplots_adjust(right=0.85)
        #cbar_ax = self.fig.add_axes([0.89, 0.6, 0.02, 0.3])
        #cb = self.fig.colorbar(im0, cax=cbar_ax)


        rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        rc('text', usetex=True)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')



        self.ax[1].set_ylabel('theta', fontsize=14)
        self.ax[2].set_ylabel('eV/A',  fontsize=14)
        self.ax[1].set_xlabel('r',     fontsize=14)
        self.ax[2].set_xlabel('r',     fontsize=14)
        n=6
        m=5
        self.ax[1].set_yticks(np.linspace(0,M,m))
        self.ax[1].set_yticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'],  fontsize=14)
        self.ax[2].set_xticks(np.linspace(0,N,n))
        self.ax[1].set_xticklabels(['%.0f'%i for i in np.linspace(0,R,n)])
        self.ax[2].set_xticklabels(['%.0f'%i for i in np.linspace(0,R,n)])
        self.ax[2].grid('on')

        plt.tight_layout()
        #plt.title(self.name)
        if self.save:
            name = "%s.png"%"".join(self.name.split())
            plt.savefig(name, format="png")
            plt.close()
            print "Saved as:", name

        else:
            if self.plotEach:
                plt.draw()
                #plt.pause(0.1)
                #plt.clf()

            if self.plotAverage:
                #plt.plot(21.5,22.5, 'g*', markersize=10)
                plt.show

#------------------------------------------------------------------------------#




if __name__=='__main__':

    singleObject = Forces('../forceFiles/m4/forces80000.txt')
    singleObject.name = 'whatever'
    singleObject.plotAverage = True
    singleObject.plotMatrix()
    plt.show()
