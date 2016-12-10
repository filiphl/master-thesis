from forceDistribution import *
from os import listdir
from matplotlib import rc
import os

folderPath = '../forceFiles/m4/'

directory = listdir(folderPath)


#print directory
r     = []
dists = {}

for f in directory:
    t = int(f[6:-4])
    print t
    N = 46
    surfN = 35
    cx=22.5
    cy=22.5
    nn=8
    bw=1
    dist = ForceDistribution(N, surfN, nn, bw, cx, cy, timeStep=t)
    dist.computeDistributions()
    dists[t]=dist
    #dist.plotDistributions()
    #plt.savefig('timeStep%06d.pdf'%t)
r = np.linspace(0,16,46)





N = 46
M = 46
R = 16

FN  = []
FN2 = []
FS  = []
FS2 = []
ds  = []

pairs = []
keys  = sorted(dists.keys())
#print keys
for i in range(len(keys)/2):
    pairs.append((keys[i], keys[-1-i]))



rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

lc = ["#5FA38E", "#3F5F7F"]

d = 0
c = 0
for pair in pairs:
    print pair
    if d < 61:
        d+=2.0
    #fig, ax = plt.subplots(1,1)
    for i in xrange(2):
        dist = dists[pair[i]]
        fn = dist.normal
        fs = dist.shear

        if i:
            FN2.append(float(sum(sum(fn))))
            FS2.append(float(sum(sum(fs))))
        else:
            FN.append(float(sum(sum(fn))))
            FS.append(float(sum(sum(fs))))

    #    output = dist.transform(fn,N,M,R,22.5,22.5)
    #    radialDist = np.mean(output,0)
    #    ax.plot(radialDist, linewidth=2, label='Time step %d'%pair[i], color=lc[i])
    #    ax.set_xticks(np.linspace(0,N,9))
    #    ax.set_xticklabels(['%.0f'%i for i in np.linspace(0,R,9)])
    #    plt.ylim([-0.01, 0.035])
    #    plt.grid('on')
    #    plt.hold('on')
#
    #plt.legend()
    #ax.set_ylabel(r"$eV/\AA$", fontsize=16)
    #ax.set_xlabel(r"$r$", fontsize=18)
    #plt.title(r'$ $Compression length: %.1f \AA'%d)
    ##if not d%3:
    #plt.savefig('compression%07d.pdf'%pair[0])

    ds.append(float(d))



print len(ds), len(FN)
fig,ax = plt.subplots(1,1)
ax.plot(ds, FN,
'-h',
color=lc[0],
linewidth=2,
markersize=7,
markerfacecolor='#82BC92',
fillstyle='full',
label='Down')

ax.hold('on')

ax.plot(ds, FN2,
'-h',
color=lc[1],
linewidth=2,
markersize=7,
markerfacecolor='#383D65',
fillstyle='full',
label='Up')

ax.legend(loc=2)
ax.grid()
ax.set_ylabel(r"Total normal force [eV/\AA]", fontsize=16)
ax.set_xlabel(r"$\Delta h$ [\AA]", fontsize=16)
plt.show()

fig,ax = plt.subplots(1,1)
ax.plot(ds, FS,  linewidth=2, color=lc[0], label='Down')
ax.hold('on')
ax.plot(ds, FS2, linewidth=2, color=lc[1], label='Up')
ax.legend(loc=2)
ax.grid()
ax.set_ylabel(r"Total shear force [eV/\AA]", fontsize=16)
ax.set_xlabel(r"$\Delta h$ [\AA]", fontsize=16)
#plt.show()






#os.system('pdftk compression*.pdf cat output CL.pdf')
#os.system('rm compression*.pdf')
plt.show()
