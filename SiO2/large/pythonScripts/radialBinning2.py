from matplotlib import pyplot as plt
import numpy as np

class smooth:

    def __init__(self, N, cm, binWidth=1, nBins=12, pBin=3):

        self.N  = N
        self.cm = cm
        self.w  = 0.5
        self.pBin = pBin
        self.nBins = nBins
        self.binWidth=binWidth
        self.f=1
        self.weights, self.bins = self.computeWeights(N, cm, self.w, self.nBins, self.pBin, self.f, self.r, self.y, self.trapzoidalIntegration)

    def y(self,x,r,cm):
        return np.sqrt(r**2-(x-cm)**2)

    def r(self,x,y,cm):
        return np.sqrt((x-cm)**2 + (y-cm)**2)


    def trapzoidalIntegration(self, F, x1, x2, N, r, cm):
        x = np.linspace(x1,x2,N)
        dx = x[1]-x[0]
        area = 0
        for i in xrange(N-1):
            area += (x[i+1]-x[i])*F( (x[i] + x[i+1])/2 , r, cm)
        return area


    def computeWeights(self, N, cm, w, nBins, pBin, f, r, y, trapzoidalIntegration):
        weights = np.zeros([N,N,nBins+1])
        bins = np.linspace(0, nBins*self.binWidth, nBins+1)
        for b in xrange(nBins):
            for i in xrange(N/2-1, N):
                for j in xrange(N/2-1, N):
                    if i>=cm:
                        if j>=cm:    # NE
                            test = False
                            if r(i-w, j-w, cm) >= bins[b] and r(i-w, j-w, cm) <= bins[b+1]:  # SW Corner in bin?
                                test = True
                                #Check integral x-limits.
                                if r(i-w, y(i-w, bins[b+1], cm)+cm, cm) <= r(i-w, j+w, cm):
                                    x1 = i-w
                                else:
                                    x1 = y(j+w, bins[b+1], cm)+cm
                                    weights[i,j,b] += x1-i+w

                                if r(i+w, y(i+w, bins[b+1], cm)+cm, cm) >= r(i+w, j-w, cm):
                                    x2 = i+w
                                else:
                                    x2 = y(j-w, bins[b+1], cm)+cm

                                if r(i+w, j+w, cm) <= bins[b+1]:    # Fully inside
                                    weights[i,j,b] += 1
                                else:
                                    weights[i,j,b] += trapzoidalIntegration(y, x1, x2, 1000, bins[b+1], cm) - (x2-x1)*(j-w-cm)
                                #print i,j,x1,x2,weights[i,j,b]
                            elif r(i+w, j+w, cm) >= bins[b] and r(i+w, j+w, cm) <= bins[b+1]: # NE Corner in bin?
                                if test:
                                    print "Already been in first!"
                                test = True
                                if r(i-w, y(i-w, bins[b], cm)+cm, cm) <= r(i-w, j+w, cm):
                                    x1 = i-w
                                else:
                                    x1 = y(j+w, bins[b], cm)+cm

                                if r(i+w, y(i+w, bins[b], cm)+cm, cm) >= r(i+w, j-w, cm):
                                    x2 = i+w
                                else:
                                    x2 = y(j-w, bins[b], cm)+cm
                                    weights[i,j,b] += (i+w-x2)
                                weights[i,j,b] += (x2-x1)*(j+w-cm) - trapzoidalIntegration(y, x1, x2, 1000, bins[b], cm)

                            elif r(i-w, j-w, cm) < bins[b] and r(i+w, j+w, cm) > bins[b+1]: #SW and NE corners are outside, but cell is inside.
                                if test:
                                    print "Already been in another!"
                                if r(i-w, y(i-w, bins[b+1], cm)+cm, cm) <= r(i-w, j+w, cm): # Cross below NW corner?
                                    x1 = i-w
                                else:
                                    x1 = y(j+w, bins[b+1], cm)+cm

                                if r(i+w, y(i+w, bins[b+1], cm)+cm, cm) >= r(i+w, j-w, cm): # Cross above SE corner?
                                    x2 = i+w
                                else:
                                    x2 = y(j-w, bins[b+1], cm)+cm
                                weights[i,j,b] -= (x2-x1)*(j+w-cm) - trapzoidalIntegration(y, x1, x2, 1000, bins[b+1], cm) + (x2<i+w)*(i+w-x2)


                                if r(i-w, y(i-w, bins[b], cm)+cm, cm) <= r(i-w, j+w, cm): # Cross below NW corner?
                                    x1 = i-w
                                else:
                                    x1 = y(j+w, bins[b], cm)+cm

                                if r(i+w, y(i+w, bins[b], cm)+cm, cm) >= r(i+w, j-w, cm): # Cross above NW corner?
                                    x2 = i+w
                                else:
                                    x2 = y(j-w, bins[b], cm)+cm
                                weights[i,j,b] -= trapzoidalIntegration(y, x1, x2, 1000, bins[b], cm) - (x2-x1)*(j-w-cm) - 1 + (x1>i-w)*(x1-(i-w))


                    # Symetry
                    weights[N-i-1, N-j-1, b] = weights[N-i-1, j, b] = weights[i, N-j-1, b] = weights[i,j,b]
        return weights, bins

    def show(self, pBin):

        fig, ax = plt.subplots()
        plt.imshow(self.weights[:,:,pBin], cmap='gray_r', interpolation='nearest', origin='lower', vmin=0, vmax=1)
        plt.colorbar()

        plt.hold('on')

        for r in self.bins[1:]:
            if r==pBin:
                circle1=plt.Circle((self.cm, self.cm), r, color='red', fill=False, linewidth=2)
            else:
                circle1=plt.Circle((self.cm, self.cm), r, color='#49848B', fill=False, linewidth=2)
            ax.add_artist(circle1)


        plt.grid('on', color='#333333', linestyle='-')
        ax.set_xticks([i+0.5 for i in xrange(-1, self.N-1)], minor=False)
        ax.set_yticks([i+0.5 for i in xrange(-1, self.N-1)], minor=False)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        plt.axis([-0.5, self.N-0.5, -0.5, self.N-0.5])
        plt.savefig('weights.pdf')
        plt.draw()

if __name__ == '__main__':

    rb = smooth(16,7.5, pBin=2)

    print rb.weights[:,:,2]
    print rb.weights.max
