from matplotlib import pyplot as plt
import numpy as np
import warnings

class smooth:

    def __init__(self, N, cx, cy, binWidth=1, nBins=12, pBin=3):

        self.N  = N
        self.cx = cx
        self.cy = cy
        self.w  = 0.5
        self.pBin = pBin
        self.nBins = nBins
        self.binWidth=binWidth
        self.f=1
        self.weights, self.bins = self.computeWeights(N, self.w, self.nBins, self.pBin, self.f, self.r, self.x, self.y, self.trapzoidalIntegration)

    def x(self,y,r):
        return np.sqrt(r**2-(y-self.cy)**2)

    def y(self,x,r):
        return np.sqrt(r**2-(x-self.cx)**2)

    def r(self,x,y):
        return np.sqrt((x-self.cx)**2 + (y-self.cy)**2)

    warnings.filterwarnings("ignore")

    def trapzoidalIntegration(self, F, x1, x2, N, r):
        x = np.linspace(x1,x2,N)
        dx = x[1]-x[0]
        area = 0
        for i in xrange(N-1):
            area += (x[i+1]-x[i])*F( (x[i] + x[i+1])/2 , r)
        return area


    def computeWeights(self, N, w, nBins, pBin, f, r, x, y, trapzoidalIntegration):
        weights = np.zeros([N,N,nBins+1])
        bins = np.linspace(0, nBins*self.binWidth, nBins+1)
        for b in xrange(nBins):
            for i in xrange(N):
                for j in xrange(N):
                    if i>=self.cx:
                        if j>=self.cy:    # NE
                            test = False
                            if r(i-w, j-w) >= bins[b] and r(i-w, j-w) <= bins[b+1]:  # SW Corner in bin?
                                test = True
                                #Check integral x-limits.
                                if r(i-w, y(i-w, bins[b+1])+self.cy) <= r(i-w, j+w):
                                    x1 = i-w
                                else:
                                    x1 = x(j+w, bins[b+1])+self.cx
                                    weights[i,j,b] += x1-i+w

                                if r(i+w, y(i+w, bins[b+1])+self.cy) >= r(i+w, j-w):
                                    x2 = i+w
                                else:
                                    x2 = x(j-w, bins[b+1])+self.cx

                                if r(i+w, j+w) <= bins[b+1]:    # Fully inside
                                    weights[i,j,b] += 1
                                else:
                                    weights[i,j,b] += trapzoidalIntegration(y, x1, x2, 1000, bins[b+1]) - (x2-x1)*(j-w-self.cy)
                                #print i,j,x1,x2,weights[i,j,b]
                            elif r(i+w, j+w) >= bins[b] and r(i+w, j+w) <= bins[b+1]: # NE Corner in bin?
                                if test:
                                    print "Already been in first!"
                                test = True
                                if r(i-w, y(i-w, bins[b])+self.cy) <= r(i-w, j+w):
                                    x1 = i-w
                                else:
                                    x1 = x(j+w, bins[b])+self.cx

                                if r(i+w, y(i+w, bins[b])+self.cy) >= r(i+w, j-w):
                                    x2 = i+w
                                else:
                                    x2 = x(j-w, bins[b])+self.cx
                                    weights[i,j,b] += (i+w-x2)
                                weights[i,j,b] += (x2-x1)*(j+w-self.cy) - trapzoidalIntegration(y, x1, x2, 1000, bins[b])

                            elif r(i-w, j-w) < bins[b] and r(i+w, j+w) > bins[b+1]: #SW and NE corners are outside, but cell is inside.
                                if test:
                                    print "Already been in another!"
                                if r(i-w, y(i-w, bins[b+1])+self.cy) <= r(i-w, j+w): # Cross below NW corner?
                                    x1 = i-w
                                else:
                                    x1 = x(j+w, bins[b+1])+self.cx

                                if r(i+w, y(i+w, bins[b+1])+self.cy) >= r(i+w, j-w): # Cross above SE corner?
                                    x2 = i+w
                                else:
                                    x2 = x(j-w, bins[b+1])+self.cx
                                weights[i,j,b] -= (x2-x1)*(j+w-self.cy) - trapzoidalIntegration(y, x1, x2, 1000, bins[b+1]) + (x2<i+w)*(i+w-x2)


                                if r(i-w, y(i-w, bins[b])+self.cy) <= r(i-w, j+w): # Cross below NW corner?
                                    x1 = i-w
                                else:
                                    x1 = x(j+w, bins[b])+self.cx

                                if r(i+w, y(i+w, bins[b])+self.cy) >= r(i+w, j-w): # Cross above NW corner?
                                    x2 = i+w
                                else:
                                    x2 = x(j-w, bins[b])+self.cx
                                weights[i,j,b] -= trapzoidalIntegration(y, x1, x2, 1000, bins[b]) - (x2-x1)*(j-w-self.cy) - 1 + (x1>i-w)*(x1-(i-w))


                    # Symetry
                    weights[2*self.cx-i, 2*self.cy-j, b] = weights[2*self.cx-i, j, b]= weights[i, 2*self.cy-j, b]= weights[i,j,b]
        return weights, bins

    def show(self, pBin):

        fig, ax = plt.subplots()
        plt.imshow(self.weights[:,:,pBin], cmap='gray_r', interpolation='nearest', origin='lower', vmin=0, vmax=1)
        plt.colorbar()

        plt.hold('on')

        for r in self.bins[1:]:
            if r==pBin:
                circle1=plt.Circle((self.cx, self.cy), r, color='red', fill=False, linewidth=2)
            else:
                circle1=plt.Circle((self.cx, self.cy), r, color='#49848B', fill=False, linewidth=2)
            ax.add_artist(circle1)


        plt.grid('on', color='#333333', linestyle='-')
        ax.set_xticks([i+0.5 for i in xrange(-1, self.N-1)], minor=False)
        ax.set_yticks([i+0.5 for i in xrange(-1, self.N-1)], minor=False)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        plt.axis([-0.5, self.N-0.5, -0.5, self.N-0.5])
        #plt.savefig('weights.pdf')
        plt.draw()


    def showInteractive(self):
        from matplotlib.widgets import Slider, Button, RadioButtons

        fig, ax = plt.subplots()
        im = plt.imshow(self.weights[:,:,0], cmap='gray_r', interpolation='nearest', origin='lower', vmin=0, vmax=1)
        plt.colorbar()

        plt.hold('on')

        for r in self.bins[1:]:
            circle1=plt.Circle((self.cy, self.cx), r, color='#49848B', fill=False, linewidth=2)
            ax.add_artist(circle1)


        plt.grid('on', color='#333333', linestyle='-')
        ax.set_xticks([i+0.5 for i in xrange(-1, self.N-1)], minor=False)
        ax.set_yticks([i+0.5 for i in xrange(-1, self.N-1)], minor=False)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        plt.axis([-0.5, self.N-0.5, -0.5, self.N-0.5])
        #plt.savefig('weights.pdf')


        #Interactive part

        axcolor = 'white'
        axBin = plt.axes([0.19, 1.5, 0.65, 0.03], axisbg=axcolor)
        sBin = Slider(axBin, 'Bin', 0, 10, valinit=0, valfmt='%0.0f', color='#49848B')


        def update(val):
            Bin = int(sBin.val)
            im.set_data(self.weights[:,:,Bin])
            fig.canvas.draw_idle()
        sBin.on_changed(update)

        upax   = plt.axes([0.445, 0.93, 0.1, 0.04])
        downax = plt.axes([0.345, 0.93, 0.1, 0.04])
        upbutton   = Button(upax,   '+', color=axcolor, hovercolor='0.9')
        downbutton = Button(downax, '-', color=axcolor, hovercolor='0.9')


        def up(event):
            if sBin.val < 8:
                sBin.val += 1
            Bin = int(sBin.val)
            im.set_data(self.weights[:,:,Bin])
            fig.canvas.draw_idle()
        upbutton.on_clicked(up)


        def down(event):
            if sBin.val > 0:
                sBin.val -= 1
            Bin = int(sBin.val)
            im.set_data(self.weights[:,:,Bin])
            fig.canvas.draw_idle()
        downbutton.on_clicked(down)

        plt.show()





if __name__ == '__main__':
    import sys

    N  = 24
    cx = 11.5
    cy = 11.5
    bw = 1.8

    if len(sys.argv)>1:
        N  = int(sys.argv[1])
        cx = float(sys.argv[2])
        cy = float(sys.argv[3])
        bw = float(sys.argv[4])
    rb = smooth(N, cx, cy, bw)
    rb.showInteractive()
