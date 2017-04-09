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
        self.radialBinning = smooth(N, cx, cy, binWidth, nBins=int(16/binWidth))
        #self.surf.plotPlanes()

    def loadSurface(self, filePath=False, N=46, s=5):
        if filePath:
            filePath = filePath.rstrip('.pkl') + 'N%ds%d.pkl'%(N,s)
            try:
                with open(filePath, 'rb') as input:
                    print "Loaded surface file with N=%d and s=%d."%(N,s)
                    return pickle.load(input)
            except:
                print "Couldn't load surface file."
        s = SurfaceRegression('../surfaceFiles/', N, False, s)
        if filePath:
            with open(filePath, 'wb') as output:
                pickle.dump(s, output, pickle.HIGHEST_PROTOCOL)
        return s

    def loadForces(self, filePath=False):
        if filePath:
            try:
                with open(filePath, 'rb') as input:
                    print "Loaded force file."
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
        Fs = np.zeros((self.N, self.N, 3))

        for i in xrange(self.N):
            for j in xrange(self.N):
                if not np.isnan( np.cos( self.surf.getAngle( self.force.matrix[i][j], self.surf.grid[ int(i/self.mapping),int(j/self.mapping) ] ) ) ):

                    a = np.dot( self.force.matrix[i,j], self.surf.grid[ int(i/self.mapping),int(j/self.mapping) ] )
                    b = np.dot(self.force.matrix[i][j], -self.surf.grid[int(i/self.mapping)][int(j/self.mapping)] )
                    if a > b:
                        self.normal[i,j] = a
                        Fs[i,j] = self.force.matrix[i,j] - self.surf.grid[ int(i/self.mapping), int(j/self.mapping) ] * self.normal[i,j]
                    else:
                        self.normal[i,j] = b
                        Fs[i,j] = self.force.matrix[i,j] + self.surf.grid[ int(i/self.mapping), int(j/self.mapping) ] * self.normal[i,j]

                    self.shear [i,j] = np.sqrt(sum(Fs[i,j]**2))
                    #x = i-self.cx
                    #y = j-self.cy
                    #self.shear[i,j] *= np.dot( self.surf.Norm(np.asarray([x,y])), Fs[i,j,:2] )
                    #angle = self.surf.getAngle( self.shear[i,j], self.surf.grid[ int(i/self.mapping), int(j/self.mapping) ] )
                    # Dot product of shear vector with unit vector from center to current point (i,j).



        self.r  = self.radialBinning.bins[1:]
        self.F  = np.zeros(self.radialBinning.nBins)
        self.FN = np.zeros(self.radialBinning.nBins)
        self.FS = np.zeros(self.radialBinning.nBins)


        for k in xrange(self.radialBinning.nBins):
            area = np.pi * (self.radialBinning.bins[k+1]**2-self.radialBinning.bins[k]**2)
            value  = 0.0
            valueN = 0.0
            valueS = 0.0
            for i in xrange(self.N):
                for j in xrange(self.N):
                    if not np.isnan(self.force.absoluteForces[i,j]):
                        value  += self.force.absoluteForces[i,j] * self.radialBinning.weights[i,j,k]
                        valueN += self.normal[i,j] * self.radialBinning.weights[i,j,k]
                        valueS += self.shear[i,j]  * self.radialBinning.weights[i,j,k]
            self.F[k]  = value /area
            self.FN[k] = valueN/area
            self.FS[k] = valueS/area


    def plotDistributions(self):
        fig, ax = plt.subplots(2,3, figsize=(15,10), sharey='row')
        fig.suptitle('Bin width = %.1f'%self.bw)

        cmax = max(self.force.absoluteForces.max(), self.normal.max(), self.shear.max())
        c=0
        for f in [self.force.absoluteForces, self.normal, self.shear]:
            im = ax[1,c].imshow(f, interpolation='nearest', cmap="hot_r", vmin=0, vmax=cmax)
            c+=1

        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.89, 0.1, 0.02, 0.3648])
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
            ax[0,c].set_xlim([0,16])
            #ax[0,c].set_ylim([0,ymax*1.05])
            ax[0,c].grid('on')

            c+=1
        ax[0,0].set_ylabel(r"$eV/\AA$", fontsize=16)
        ax[0,0].set_title(r"$ $Absolute", fontsize=16)
        ax[0,1].set_title(r"$ $Normal", fontsize=16)
        ax[0,2].set_title(r"$ $Shear", fontsize=16)




#------------------------------------------------------------------------------#

if __name__ == '__main__':

    for bw in np.linspace(1,1,1):
        N = 46
        surfN = 30
        #bw = 1.2
        cx=22.5
        cy=22.5
        nn=8
        dist = ForceDistribution(N, surfN, nn, bw, cx, cy)
        dist.computeDistributions()
        dist.plotDistributions()


    plt.show()
