from surface import *
from compareMatrices import *
from radialBinning import smooth
from matplotlib import rc
from scipy import ndimage
import matplotlib as mpl
import pickle

class ForceDistribution:

    def __init__(self, N, surfN, nearestNeighbor, binWidth, cx, cy, timeStep=False):
        self.N = N
        self.surfN = surfN
        self.nn = nearestNeighbor
        self.bw = binWidth
        self.cx=cx
        self.cy=cy
        self.mapping = float(N)/surfN

        self.timeStep = timeStep    # Only analyse single frame
        self.surf  = self.loadSurface('../dataFiles/surface.pkl', N=self.surfN, s=self.nn)
        self.force = self.loadForces('../dataFiles/forces.pkl')
        #self.radialBinning = smooth(N, cx, cy, binWidth, nBins=int(16/binWidth))
        #self.surf.plotPlanes()

    def loadSurface(self, filePath=False, N=46, s=5):

        if filePath:
            if self.timeStep:
                filePath = filePath.rstrip('.pkl') + '_t%d_N%ds%d.pkl'%(self.timeStep,N,s)
            else:
                filePath = filePath.rstrip('.pkl') + '_N%ds%d.pkl'%(N,s)

            try:
                with open(filePath, 'rb') as input:
                    print "Loaded surface file", filePath
                    pkl = pickle.load(input)
                    return pkl
            except:
                print "Couldn't load surface file."

        if self.timeStep:
            s = SurfaceRegression('../surfaceFiles/Surface%d/'%self.timeStep, N, False, s)
        else:
            s = SurfaceRegression('../surfaceFiles/', N, False, s)

        if filePath:
            with open(filePath, 'wb') as output:
                pickle.dump(s, output, pickle.HIGHEST_PROTOCOL)

        return s

    def loadForces(self, filePath=False):
        if filePath:
            if self.timeStep:
                filePath = filePath.rstrip('.pkl') + '_t%d.pkl'%self.timeStep
            try:
                with open(filePath, 'rb') as input:
                    pkl = pickle.load(input)
                    print "Loaded force file  ", filePath
                    return pkl
            except:
                print "Couldn't load force file."

        if self.timeStep:
            F = Forces('../forceFiles/partitions/forces%d.txt'%self.timeStep)
        else:
            F = Forces('../forceFiles/forcesAll.txt')

        F.plotAverage = True
        F.name = 'Averaged normal force'

        if filePath:
            with open(filePath, 'wb') as output:
                pickle.dump(F, output, pickle.HIGHEST_PROTOCOL)

        return F

    def transform(self, matrix, N, M, R, cx, cy):

        r = np.linspace(0,R,N)
        theta = np.linspace(0, 2*np.pi, M)
        myevalmatrix = np.zeros((N, M, 2))
        for i in range(N):
            for j in range(M):
                myevalmatrix[i, j,:] = np.asarray([cx+r[i]*np.cos(theta[j]), cy+r[i]*np.sin(theta[j])])

        return ndimage.map_coordinates(matrix, np.transpose(myevalmatrix[:, :]), order=1)

#------------------------------------------------------------------------------#
    def computeDistributions(self):


    #surf.plotPlanes()
    #plt.figure()
    #force.plotMatrix()
    #plt.show()
    #radialBinning.show(3)

        self.normal = copy.deepcopy(self.force.absoluteForces)
        self.shear  = copy.deepcopy(self.force.absoluteForces)
        Fs = np.zeros((self.N, self.N, 3))

        for i in xrange(self.N):
            for j in xrange(self.N):
                if not np.isnan( np.cos( self.surf.getAngle( self.force.matrix[i][j], self.surf.grid[ int(i/self.mapping),int(j/self.mapping) ] ) ) ):

                    a = np.dot( self.force.matrix[i,j],  self.surf.grid[ int(i/self.mapping),int(j/self.mapping) ] )
                    b = np.dot( self.force.matrix[i,j], -self.surf.grid[ int(i/self.mapping),int(j/self.mapping) ] )

                    if a > b:
                        self.normal[i,j] = b # a
                        Fs[i,j] = self.force.matrix[i,j] - self.surf.grid[ int(i/self.mapping), int(j/self.mapping) ] * self.normal[i,j]
                    else:
                        self.normal[i,j] = b
                        Fs[i,j] = self.force.matrix[i,j] + self.surf.grid[ int(i/self.mapping), int(j/self.mapping) ] * self.normal[i,j]

                    self.shear [i,j] = np.sqrt(sum(Fs[i,j]**2))
                    x = i-self.cx
                    y = j-self.cy
                    self.shear[i,j] *= np.dot( self.surf.Norm(np.asarray([x,y])), self.surf.Norm(Fs[i,j,:2]) )
                    #angle = self.surf.getAngle( self.shear[i,j], self.surf.grid[ int(i/self.mapping), int(j/self.mapping) ] )
                    # Dot product of shear vector with unit vector from center to current point (i,j).





    def plotDistributions(self):
        fig, ax = plt.subplots(2,3, figsize=(15,10), sharey='row', sharex='col')
        cmax = max(self.force.absoluteForces.max(), self.normal.max(), self.shear.max())

        N = 46
        M = 46
        R = 16


        c=0
        for f in [self.force.absoluteForces, self.normal, self.shear]:
            output = self.transform(f,N,M,R,22.5,22.5)
            im=ax[0,c].pcolor(output, vmin=0, vmax=cmax)

            radialDist = np.mean(output,0)
            ax[1,c].plot(radialDist, linewidth=2, color="#478684")
            ax[0,c].set_xlim([0,46])
            ax[1,c].set_xlabel(r"$r$", fontsize=18)
            ax[1,c].set_xticks(np.linspace(0,N,6))
            ax[1,c].set_xticklabels(['%.0f'%i for i in np.linspace(0,R,6)])
            #ax[0,c].set_ylim([0,ymax*1.05])
            ax[1,c].grid('on')
            ax[1,c].set_ylim([-0.01,0.025])
            c+=1

        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.89, 0.535, 0.02, 0.3648])
        fig.colorbar(im, cax=cbar_ax, format=ticker.FuncFormatter(self.force.fmt))


        rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        rc('text', usetex=True)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')


        ax[0,0].set_ylim([0,46])
        ax[0,0].set_ylabel(r"$\theta$", fontsize=18)
        ax[1,0].set_ylabel(r"$eV/\AA$", fontsize=16)
        ax[0,0].set_title(r"$ $Magnitude", fontsize=16)
        ax[0,1].set_title(r"$ $Normal", fontsize=16)
        ax[0,2].set_title(r"$ $Shear", fontsize=16)

        ax[0,0].set_yticks(np.linspace(0,M,5))
        ax[0,0].set_yticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'],  fontsize=14)

        if self.timeStep:
            fig.suptitle(r"$ $Time step %d"%self.timeStep, fontsize=16)

#------------------------------------------------------------------------------#

if __name__ == '__main__':

    for bw in np.linspace(1,1,1):
        N = 46
        surfN = 23
        #bw = 1.2
        cx=22.5
        cy=22.5
        nn=5
        dist = ForceDistribution(N, surfN, nn, bw, cx, cy, timeStep=False)
        dist.computeDistributions()
        dist.plotDistributions()


    plt.show()
