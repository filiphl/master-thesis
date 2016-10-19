from matplotlib import pyplot as plt
import numpy as np

'''
class smooth:

    def __init__(self,N,cm):

        self.N  = N
        self.cm = cm
        self.w  = 0.5

        matrix = np.zeros([N,N])

    def y(x,r,cm):
        return np.sqrt(r**2-(x-cm)**2)

    def R(x,y,cm):
        return np.sqrt((x-cm)**2 + (y-cm)**2)





'''

def trapzoidalIntegration(F,x1,x2,N, r, cm):
    x = np.linspace(x1,x2,N)
    dx = x[1]-x[0]
    area = 0
    for i in xrange(N-1):
        area += (x[i+1]-x[i])*F( (x[i] + x[i+1])/2 , r, cm)
    return area

N  = 10
cm = 4.5
w  = 0.5

matrix = np.zeros([N,N])
inside = []

def y(x,r,cm):
    return np.sqrt(r**2-(x-cm)**2)

def r(x,y,cm):
    return np.sqrt((x-cm)**2 + (y-cm)**2)

nBins = 10
weights = np.zeros([N,N,nBins+1])
bins = np.linspace(0,nBins,nBins+1)
for b in xrange(nBins):
    for i in xrange(N/2-1, N):
        for j in xrange(N/2-1, N):
            if i>cm:
                if j>cm:    # NE
                    if r(i-w, j-w, cm) >= bins[b]:  # SW Corner in bin?
                        if r(i-w, j-w, cm) < bins[b+1]:
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
                            weights[i,j,b] += trapzoidalIntegration(y, x1, x2, 100, bins[b+1], cm) - (x2-x1)*(j-w-cm)
                            #print i,j,x1,x2,weights[i,j,b]
                    if r(i+w, j+w, cm) >= bins[b]:  # NE Corner in bin?
                        if r(i+w, j+w, cm) < bins[b+1]:
                            if r(i-w, y(i-w, bins[b], cm)+cm, cm) <= r(i-w, j+w, cm):
                                x1 = i-w
                            else:
                                x1 = y(j+w, bins[b], cm)+cm

                            if r(i+w, y(i+w, bins[b], cm)+cm, cm) >= r(i+w, j-w, cm):
                                x2 = i+w
                            else:
                                x2 = y(j-w, bins[b], cm)+cm
                                weights[i,j,b] += (i+w-x2)
                            weights[i,j,b] += (x2-x1)*(j+w-cm) - trapzoidalIntegration(y, x1, x2, 100, bins[b], cm)
            weights[N-i-1, N-j-1, b] = weights[N-i-1, j, b] = weights[i, N-j-1, b] = weights[i,j,b]



        '''
                if r(i+w, j+w, cm) >= 2:
                    if r(i+w, j+w, cm) < 3:
                        matrix[i,j] = 1
                        inside.append([i,j])

            else:   #j<=cm
                if r(i-w, j+w, cm) >= 2:
                    if r(i-w, j+w, cm) < 3:
                        matrix[i,j] = 1
                        inside.append([i,j])
                if r(i+w, j-w, cm) >= 2:
                    if r(i+w, j-w, cm) < 3:
                        matrix[i,j] = 1
                        inside.append([i,j])

        else:       # i<cm
            if j<cm:
                if r(i+w, j+w, cm) >= 2:
                    if r(i+w, j+w, cm) < 3:
                        matrix[i,j] = 1
                        inside.append([i,j])
                if r(i-w, j-w, cm) >= 2:
                    if r(i-w, j-w, cm) < 3:
                        matrix[i,j] = 1
                        inside.append([i,j])

            else:
                if r(i+w, j-w, cm) >= 2:
                    if r(i+w, j-w, cm) < 3:
                        matrix[i,j] = 1
                        inside.append([i,j])
                if r(i-w, j+w, cm) >= 2:
                    if r(i-w, j-w, cm) < 3:
                        matrix[i,j] = 1
                        inside.append([i,j])
        '''
c=1
for i in inside:
    matrix[i[0],i[1]]=c
    c+=1



fig, ax = plt.subplots()
plt.imshow(weights[:,:,3], cmap='gray_r', interpolation='nearest', origin='lower', vmin=0, vmax=1)
plt.colorbar()

plt.hold('on')

for r in range(1,7):
    circle1=plt.Circle((cm, cm), r, color='#49848B', fill=False, linewidth=2)
    ax.add_artist(circle1)


plt.grid('on', color='#333333', linestyle='-')
ax.set_xticks([i+0.5 for i in xrange(-1,N-1)], minor=False)
ax.set_yticks([i+0.5 for i in xrange(-1,N-1)], minor=False)
ax.set_xticklabels([])
ax.set_yticklabels([])
plt.axis([-0.5, N-0.5, -0.5, N-0.5])
#plt.savefig('weights.pdf')
plt.show()
