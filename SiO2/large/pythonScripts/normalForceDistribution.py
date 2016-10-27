from surface import *
from compareMatrices import *
from radialBinning import smooth
from matplotlib import rc
import pickle

class ForceDistribution:

    def __init__(self, N, surfN, nearestNeighbor, binWidth, cx, cy):
        self.N = N
        self.surfN = surfN
        self.nn = nearestNeighbor
        self.bw = binWidth
        self.cx=cx
        self.cy=cy
        self.mapping = float(N)/surfN

        self.surf  = self.loadSurface('../dataFiles/surface.pkl', N=surfN, s=nn)
        self.force = self.loadForces('../dataFiles/forces.pkl')
        self.radialBinning = smooth(N, cx, cy, binWidth, nBins=int(12/binWidth))

    def loadSurface(self, filePath=False, N=45, s=5):
        if filePath:
            filePath = filePath.rstrip('.pkl') + 'N%ds%d.pkl'%(N,s)
            try:
                with open(filePath, 'rb') as input:
                    print "Loaded surface file with N=%d and s=%d."%(N,s)
                    return pickle.load(input)
            except:
                print "Couldn't load surface file."
        s = SurfaceRegression('../dumpFiles/', N, False, s)
        if filePath:
            with open(filePath, 'wb') as output:
                pickle.dump(s, output, pickle.HIGHEST_PROTOCOL)
        return s
    #surf.plotPlanes()

    def loadForces(self, filePath=False):
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
    def computeDistributions(self):


    #surf.plotPlanes()
    #plt.figure()
    #force.plotMatrix()
    #plt.show()
    #radialBinning.show(3)

        self.normal = copy.deepcopy(self.force.absoluteForces)
        self.shear  = copy.deepcopy(self.force.absoluteForces)
        for i in xrange(45):
            for j in xrange(45):
                if not np.isnan( np.cos( self.surf.getAngle( self.force.matrix[i][j], self.surf.grid[int(i/self.mapping)][int(j/self.mapping)] ) ) ):
                    angle = self.surf.getAngle( self.force.matrix[i][j], self.surf.grid[int(i/self.mapping)][int(j/self.mapping)] )
                    self.normal[i][j] *= np.cos(angle)
                    self.shear [i][j] *= np.sin(angle)


        self.r  = self.radialBinning.bins[1:]
        self.F  = np.zeros(self.radialBinning.nBins)
        self.FN = np.zeros(self.radialBinning.nBins)
        self.FS = np.zeros(self.radialBinning.nBins)


        for k in xrange(self.radialBinning.nBins):
            area = np.pi * (self.radialBinning.bins[k+1]**2-self.radialBinning.bins[k]**2)
            value  = 0.0
            valueN = 0.0
            valueS = 0.0
            for i in xrange(45):
                for j in xrange(45):
                    if not np.isnan(self.force.absoluteForces[i,j]):
                        value  += self.force.absoluteForces[i,j] * self.radialBinning.weights[i,j,k]
                        valueN += self.normal[i,j] * self.radialBinning.weights[i,j,k]
                        valueS += self.shear[i,j]  * self.radialBinning.weights[i,j,k]
            self.F[k]  = value /area
            self.FN[k] = valueN/area
            self.FS[k] = valueS/area


    def plotDistributions(self):
        fig, ax = plt.subplots(2,3, figsize=(15,10))
        fig.suptitle('Bin width = %.1f'%self.bw)

        cmax = max(self.force.absoluteForces.max(), self.normal.max(), self.shear.max())
        c=0
        for f in [self.force.absoluteForces, self.normal, self.shear]:
            im = ax[1,c].imshow(f, interpolation='nearest', cmap="hot_r", vmin=0, vmax=cmax)
            c+=1

        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.9, 0.1, 0.02, 0.35])
        fig.colorbar(im, cax=cbar_ax, format=ticker.FuncFormatter(self.force.fmt))


# ------------------------ Radial distributions ------------------------------ #


        rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        rc('text', usetex=True)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        ymax = max(self.F.max(), self.FN.max(), self.FS.max())
        c=0
        for f in [self.F, self.FN, self.FS]:

            ax[0,c].plot(self.r,f,
             '-h',
             color="#478684",
             linewidth=2,
             markersize=7,
             markerfacecolor='#3b617c',
             fillstyle='full')
            ax[0,c].set_xlabel(r"$r$", fontsize=18)
            ax[0,c].set_xlim([0,12])
            ax[0,c].set_ylim([0,ymax*1.05])
            ax[0,c].grid('on')

            c+=1
        ax[0,0].set_ylabel(r"$eV/\AA$", fontsize=16)
        ax[0,0].set_title(r"$ $Absolute", fontsize=16)
        ax[0,1].set_title(r"$ $Normal", fontsize=16)
        ax[0,2].set_title(r"$ $Shear", fontsize=16)




#------------------------------------------------------------------------------#

if __name__ == '__main__':

    for bw in np.linspace(1,1,1):
        N = 45
        surfN = 11
        nn = 0
        #bw = 1.2
        cx=22
        cy=22
        nn=5
        dist = ForceDistribution(N, surfN, nn, bw, cx, cy)
        dist.computeDistributions()
        dist.plotDistributions()

    plt.show()
