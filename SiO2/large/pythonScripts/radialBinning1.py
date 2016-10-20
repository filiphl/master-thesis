from matplotlib import pyplot as plt
import numpy as np

class smooth:

    def __init__(self,N,cm):

        self.N  = N
        self.cm = cm
        self.w  = 0.5
        self.pBin = 10
        self.nBins = N
        f=1
        self.weights = np.zeros([N,N,self.nBins+1])
        self.bins = np.linspace(0,self.nBins, f*self.nBins+1)

        self.computeWeights(N, cm, self.w, self.nBins, self.pBin, self.bins)




    def trapzoidalIntegration(self, F,x1,x2,N, r):
        x = np.linspace(x1,x2,N)
        dx = x[1]-x[0]
        area = 0
        for i in xrange(N-1):
            area += (x[i+1]-x[i])*F( (x[i] + x[i+1])/2 , r)
        return area


    def y(self, x,r):
        return np.sqrt(r**2-(x-self.cm)**2)

    def r(self, x,y):
        return np.sqrt((x-self.cm)**2 + (y-self.cm)**2)





    def computeWeights(self, N, cm, w, nBins, pBin, bins):

        for b in xrange(nBins):
            for i in xrange(N/2-1, N):
                for j in xrange(N/2-1, N):
                    if i>cm:
                        if j>cm:    # NE
                            if self.r(i-w, j-w) >= bins[b]:  # SW Corner in bin?
                                if self.r(i-w, j-w) < bins[b+1]:
                                    #Check integral x-limits.
                                    if self.r(i-w, self.y(i-w, bins[b+1])+cm) <= self.r(i-w, j+w):
                                        x1 = i-w
                                    else:
                                        x1 = self.y(j+w, bins[b+1])+cm
                                        self.weights[i,j,b] += x1-i+w

                                    if self.r(i+w, self.y(i+w, bins[b+1])+cm) >= self.r(i+w, j-w):
                                        x2 = i+w
                                    else:
                                        x2 = self.y(j-w, bins[b+1])+cm
                                    self.weights[i,j,b] += self.trapzoidalIntegration(self.y, x1, x2, 100, bins[b+1]) - (x2-x1)*(j-w-cm)
                                    #print i,j,x1,x2,self.weights[i,j,b]
                            if self.r(i+w, j+w) >= bins[b]:  # NE Corner in bin?
                                if self.r(i+w, j+w) < bins[b+1]:
                                    if self.r(i-w, self.y(i-w, bins[b])+cm) <= self.r(i-w, j+w):
                                        x1 = i-w
                                    else:
                                        x1 = self.y(j+w, bins[b])+cm

                                    if self.r(i+w, self.y(i+w, bins[b])+cm) >= self.r(i+w, j-w):
                                        x2 = i+w
                                    else:
                                        x2 = self.y(j-w, bins[b])+cm
                                        self.weights[i,j,b] += (i+w-x2)
                                    self.weights[i,j,b] += (x2-x1)*(j+w-cm) - self.trapzoidalIntegration(self.y, x1, x2, 100, bins[b])

                    if self.r(i-w, j-w) < bins[b] and self.r(i+w, j+w) > bins[b+1]:

                        if self.r(i-w, self.y(i-w, bins[b+1])+cm) <= self.r(i-w, j+w):
                            x1 = i-w
                        else:
                            x1 = self.y(j+w, bins[b+1])+cm

                        if self.r(i+w, self.y(i+w, bins[b+1])+cm) >= self.r(i+w, j-w):
                            x2 = i+w
                        else:
                            x2 = self.y(j-w, bins[b+1])+cm
                        self.weights[i,j,b] -= (x2-x1)*(j+w-cm) - self.trapzoidalIntegration(self.y, x1, x2, 1000, bins[b+1]) + (x2<i+w)*(i+w-x2)


                        if self.r(i-w, self.y(i-w, bins[b])+cm) <= self.r(i-w, j+w):
                            x1 = i-w
                        else:
                            x1 = self.y(j+w, bins[b])+cm

                        if self.r(i+w, self.y(i+w, bins[b])+cm) >= self.r(i+w, j-w):
                            x2 = i+w
                        else:
                            x2 = self.y(j-w, bins[b])+cm
                        self.weights[i,j,b] -= self.trapzoidalIntegration(self.y, x1, x2, 1000, bins[b]) - (x2-x1)*(j-w-cm) - 1 + (x1>i-w)*(x1-(i-w))
                        if b == pBin:
                            print i,j,self.weights[i,j,b]

                    # Symetry
                    self.weights[N-i-1, N-j-1, b] = self.weights[N-i-1, j, b] = self.weights[i, N-j-1, b] = self.weights[i,j,b]






if __name__ == '__main__':
    re = smooth(10,4.5)



    fig, ax = plt.subplots()
    plt.imshow(re.weights[:,:,re.pBin], colormap='gray_r', interpolation='nearest', origin='lower', vmin=0, vmax=1)
    plt.colorbar()

    plt.hold('on')

    for r in re.bins[1:]:
        if r==re.pBin:
            circle1=plt.Circle((re.cm), r, color='red', fill=False, linewidth=2)
        else:
            circle1=plt.Circle((re.cm), r, color='#49848B', fill=False, linewidth=2)
        ax.add_artist(circle1)


    plt.grid('on', color='#333333', linestyle='-')
    ax.set_xticks([i+0.5 for i in xrange(-1,re.N-1)], minor=False)
    ax.set_yticks([i+0.5 for i in xrange(-1,re.N-1)], minor=False)
    #ax.set_xticklabels([])
    #ax.set_yticklabels([])
    plt.axis([-0.5, re.N-0.5, -0.5, re.N-0.5])
    #plt.savefig('weights.pdf')
    plt.show()
