from forceDistribution import *
from os import listdir
from matplotlib import rc


folderPath = '../forceFiles/partitions10000/'

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

pairs = []
keys = sorted(dists.keys())[:-1]

print keys
for i in range(len(keys)/2):
    pairs.append((keys[i], keys[-1-i]))



rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

lc = ['#5fa38e', '#3f5f7f']
d = 0
for pair in pairs:
    d+=5
    fig, ax = plt.subplots(1,1)
    for i in xrange(2):
        dist = dists[pair[i]]
        f = dist.normal
        output = dist.transform(f,N,M,R,22.5,22.5)
        radialDist = np.mean(output,0)
        ax.plot(radialDist, linewidth=2, label='Time step %d'%pair[i], color=lc[i])
        ax.set_xticks(np.linspace(0,N,9))
        ax.set_xticklabels(['%.0f'%i for i in np.linspace(0,R,9)])
        plt.ylim([-0.005, 0.025])
        plt.grid('on')
        plt.hold('on')
    plt.legend()

    ax.set_ylabel(r"$eV/\AA$", fontsize=16)
    ax.set_xlabel(r"$r$", fontsize=18)


    plt.title(r'$ $Compression length: %d \AA'%d)
    plt.savefig('compressionLength%02d.pdf'%d)
    plt.show()
