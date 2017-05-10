from surface import *
from compareMatrices import *
from radialBinning import smooth
from matplotlib import rc
from scipy import ndimage
import matplotlib as mpl
import pickle

class ForceDistribution:

    def __init__(self, N, surfN, nearestNeighbor, binWidth, cx, cy, timeStep=False, verbose=False):
        self.N = N
        self.surfN = surfN
        self.nn = nearestNeighbor
        self.bw = binWidth
        self.cx=cx
        self.cy=cy
        self.mapping = float(N)/surfN

        self.timeStep = timeStep    # Only analyse single frame
        self.surf  = self.loadSurface('../dataFiles/m4/surface.pkl', N=self.surfN, s=self.nn, verbose=verbose)
        self.force = self.loadForces('../dataFiles/m4/forces.pkl', verbose=verbose)

        #self.radialBinning = smooth(N, cx, cy, binWidth, nBins=int(16/binWidth))
        #self.surf.plotPlanes()

    def loadSurface(self, filePath=False, N=46, s=5, verbose=False):
        if filePath:
            if self.timeStep != 0:
                filePath = filePath.rstrip('.pkl') + '_t%d_N%ds%d.pkl'%(self.timeStep,N,s)
            else:
                filePath = filePath.rstrip('.pkl') + '_N%ds%d.pkl'%(N,s)

            try:
                with open(filePath, 'rb') as input:
                    pkl = pickle.load(input)
                    if verbose:
                        print "Loaded surface file", filePath
                    return pkl
            except:
                if verbose:
                    print "Couldn't load surface file."

        if self.timeStep:
            s = SurfaceRegression('../surfaceFiles/m4/Surface%d/'%self.timeStep, N, False, s)
        else:
            s = SurfaceRegression('../surfaceFiles/m4/', N, False, s)

        if filePath:
            with open(filePath, 'wb') as output:
                pickle.dump(s, output, pickle.HIGHEST_PROTOCOL)

        return s

    def loadForces(self, filePath=False, verbose=False):
        if filePath:
            if self.timeStep:
                filePath = filePath.rstrip('.pkl') + '_t%d.pkl'%self.timeStep
            try:
                with open(filePath, 'rb') as input:
                    pkl = pickle.load(input)
                    if verbose:
                        print "Loaded force file  ", filePath
                    return pkl
            except:
                if verbose:
                    print "Couldn't load force file."

        if self.timeStep:
            F = Forces('../forceFiles/m4/forces%d.txt'%self.timeStep, self.cx, self.cy)
        else:
            F = Forces('../forceFiles/m4/forcesAll.txt', self.cx, self.cy)

        F.plotAverage = True
        F.name = 'Averaged normal force'

        if filePath:
            with open(filePath, 'wb') as output:
                pickle.dump(F, output, pickle.HIGHEST_PROTOCOL)

        return F


    def transform(self, matrix, N, M, R):

        r = np.linspace(0,R,N)
        theta = np.linspace(0, 2*np.pi, M)
        myevalmatrix = np.zeros((N, M, 2))
        for i in range(N):
            for j in range(M):
                myevalmatrix[i, j,:] = np.asarray([self.cx+r[i]*np.cos(theta[j]), self.cy+r[i]*np.sin(theta[j])])

        return ndimage.map_coordinates(matrix, np.transpose(myevalmatrix[:, :]), order=1)

#------------------------------------------------------------------------------#
    def computeDistributions(self):
        '''The if-test ensures a positive component for the normal force, which
        is incorrect! This must be figured out. Also, the definition of shear
        force must be applied.'''

        self.normal = np.zeros(np.shape(self.force.absoluteForces))
        self.shear  = np.zeros(np.shape(self.force.absoluteForces))
        Fs = np.zeros((self.N, self.N, 3))

        for i in xrange(self.N):
            for j in xrange(self.N):
                mi = int(i/self.mapping)
                mj = int(j/self.mapping)
                self.surf.grid[mi,mj] /= np.linalg.norm(self.surf.grid[mi,mj])

                if not np.isnan( np.cos( self.surf.getAngle( self.force.matrix[i][j], self.surf.grid[mi, mj] ) ) ):
                    a = np.dot( self.force.matrix[i,j],  self.surf.grid[mi, mj] )
                    b = np.dot( self.force.matrix[i,j], -self.surf.grid[mi, mj] )
                    if a > b:
                        self.normal[i,j] = b # a
                        #self.shear[i,j] = np.sqrt(self.force.absoluteForces[i,j]**2-b**2)
                        #Fs[i,j] = self.force.matrix[i,j] + self.surf.grid[ int(i/self.mapping), int(j/self.mapping) ] * b
                    else:
                        self.normal[i,j] = b #b
                        #self.shear[i,j] = np.sqrt(self.force.absoluteForces[i,j]**2-b**2) #a

                    Fn = self.surf.grid[mi, mj]*self.normal[i,j]
                    Fs[i,j] = self.force.matrix[i,j] + Fn # The reason for the addition is that Fn points downwards, rather than upwards. Too late to change now.

                    self.shear [i,j] = np.sqrt(sum(Fs[i,j]**2))
                    x = i-self.cx
                    y = j-self.cy
                    self.shear[i,j] *= np.dot( self.surf.Norm(np.asarray([x,y])), self.surf.Norm(Fs[i,j,:2]) )
                    ##angle = self.surf.getAngle( self.shear[i,j], self.surf.grid[ int(i/self.mapping), int(j/self.mapping) ] )
                    ## Dot product of shear vector with unit vector from center to current point (i,j).





    def plotDistributions(self, onlyBottom):


        majorFontSize = 30
        minorFontSize = 26
        ticksfontSize = 26


        rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        rc('text', usetex=True)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')


        if onlyBottom:
            fig, ax = plt.subplots(1,3, figsize=(15,5.4), sharey='row')
            cmax = max(self.force.absoluteForces.max(), self.normal.max(), self.shear.max())
            N = 46
            M = 46
            R = 20


            myMap = sns.cubehelix_palette(80, start=.5, rot=-.75)
            cm = mpl.colors.ListedColormap(myMap)
            c=0
            for f in [self.force.absoluteForces, self.normal, self.shear]:
                output = self.transform(f,N,M,R)
                radialDist = np.mean(output,0)
                ax[c].plot(radialDist, linewidth=3, color="#478684")
                ax[c].set_xlabel(r'{\fontsize{36pt}{3em}\selectfont{}{$r$}' + r'{\fontsize{26pt}{3em}\selectfont{} $[7.12\mathring{A}]$}', labelpad=20)
                ax[c].set_xticks(np.linspace(0,N,6))
                ax[c].set_xticklabels(['' for i in range(6)])
                ax[c].set_xticklabels(['$%.0f$'%i for i in np.linspace(0,R,6)], fontsize=ticksfontSize)
                ax[c].grid('on')
                ax[c].set_ylim([-0.01,0.05])
                c+=1




            ax[0].set_ylabel(r"$eV/\mathring{A}$", fontsize=majorFontSize)

            yax = [r'$%.02f$'%i for i in np.arange(-0.01, 0.06, 0.01)]
            xax = [r'$%d$'%i for i in range(0,21,4)]

            ax[0].set_yticklabels(yax, fontsize=ticksfontSize)


            plt.gcf().subplots_adjust(bottom=0.23)
            #plt.tight_layout()
            if self.timeStep:
                #fig.suptitle(r"$ $Time step %d"%self.timeStep, fontsize=16)
                plt.savefig('timeSteps/radialOnly/timestep%06d_bottom.pdf'%self.timeStep)








        else:
            fig, ax = plt.subplots(2,3, figsize=(15,10), sharey='row', sharex='col')
            cmax = max(self.force.absoluteForces.max(), self.normal.max(), self.shear.max())
            N = 46
            M = 46
            R = 20


            myMap = sns.cubehelix_palette(80, start=.5, rot=-.75)
            cm = mpl.colors.ListedColormap(myMap)
            c=0
            for f in [self.force.absoluteForces, self.normal, self.shear]:
                output = self.transform(f,N,M,R)
                im=ax[0,c].pcolor(output, vmin=-0.01, vmax=0.07, cmap=cm)

                radialDist = np.mean(output,0)
                ax[1,c].plot(radialDist, linewidth=3, color="#478684")
                ax[0,c].set_xlim([0,46])
                ax[1,c].set_xlabel(r"$r$", fontsize=majorFontSize+2)
                ax[1,c].set_xticks(np.linspace(0,N,6))
                ax[1,c].set_xticklabels(['$%.0f$'%i for i in np.linspace(0,R,6)], fontsize=ticksfontSize)
                #ax[0,c].set_ylim([0,ymax*1.05])
                ax[1,c].grid('on')
                ax[1,c].set_ylim([-0.01,0.05])
                c+=1


            fig.subplots_adjust(right=0.85)
            cbar_ax = fig.add_axes([0.89, 0.535, 0.02, 0.3648])
            cbar = fig.colorbar(im, cax=cbar_ax, format='$%.02f$')#, format=ticker.FuncFormatter(self.force.fmt))
            cbar.ax.tick_params(labelsize=ticksfontSize)
            cbar.set_label(r'$eV/\mathring{A}$',size=majorFontSize, labelpad=10)

            rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
            rc('text', usetex=True)
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')


            ax[0,0].set_ylim([0,46])
            ax[0,0].set_ylabel(r"$\theta$", fontsize=majorFontSize+2)
            ax[1,0].set_ylabel(r"$eV/\mathring{A}$", fontsize=majorFontSize)
            ax[0,0].set_title(r"$ $Magnitude", fontsize=ticksfontSize)
            ax[0,1].set_title(r"$ $Normal", fontsize=ticksfontSize)
            ax[0,2].set_title(r"$ $Shear", fontsize=ticksfontSize)

            yax = [r'$%.02f$'%i for i in np.arange(-0.01, 0.06, 0.01)]
            xax = [r'$%d$'%i for i in range(0,21,4)]

            ax[0,0].set_yticks(np.linspace(0,M,5))
            ax[0,0].set_yticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'],  fontsize=ticksfontSize)
            ax[1,0].set_yticklabels(yax, fontsize=ticksfontSize)

            if self.timeStep:
                #fig.suptitle(r"$ $Time step %d"%self.timeStep, fontsize=16)
                plt.savefig('timeSteps/timestep%06d.pdf'%self.timeStep)
    #------------------------------------------------------------------------------#

if __name__ == '__main__':
    import seaborn as sns
    import sys
    for bw in np.linspace(1,1,1):
        N = 46
        surfN = 46
        #bw = 1.2
        cx=22.5
        cy=22.5
        nn=8
        a = 120000; b = 125000; c=5000
        for i in xrange(a, b, c):
            print '\r %.0f%% complete'%((i-a)*100.0/(b-a)),
            sys.stdout.flush()
            dist = ForceDistribution(N, surfN, nn, bw, cx, cy, timeStep=i, verbose=False)
            dist.computeDistributions()
            dist.plotDistributions(onlyBottom=True)

    plt.show()
