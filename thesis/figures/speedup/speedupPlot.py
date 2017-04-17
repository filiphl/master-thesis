from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np
from sys import argv



p = [2**i for i in range(7)]
#t = [61.55, 32.53, 18.16, 10.35, 6.09, 3.33, 2.02]
t = [3741, 1883, 1061, 634, 378.5, 208.25, 120.5]

speedup = np.zeros(7)

for i in xrange(len(t)):
    t[i] = int(t[i])*60 + t[i] - int(t[i])
    speedup[i] = float(t[0])/t[i]


fig, ax = plt.subplots()
'''
ax.set_xticks([0]+[2**i for i in range(1,7)])
ax.set_yticks([0]+[2**i for i in range(1,7)])
ax.set_xticklabels(['0', '2', '4', '8', '16', '32', '64'])
ax.set_yticklabels(['0', '2', '4', '8', '16', '32', '64'])
plt.axis([0, 64, 0, 64])
'''


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
ax.set_xticks([2**i for i in range(0,7)])
ax.set_yticks([2**i for i in range(0,7)])
ax.set_xticklabels([r'$%d$'%2**i for i in range(0,7)], fontsize=20)
ax.set_yticklabels([r'$%d$'%2**i for i in range(0,7)], fontsize=20)

plt.plot(p,p,
 '--',
 color="#DB7477",
 linewidth=3,
 )

plt.hold('on')

plt.plot(p,speedup,
 '-h',
 color="#478684",
 linewidth=4,
 markersize=9,
 markerfacecolor='#3b617c',
 fillstyle='full')






plt.axis([0, 64, 0, 64])
plt.xlabel("Number of processors", fontsize=20)
plt.ylabel(r"Speedup", fontsize=20)
plt.grid('on')




ax.xaxis.set_label_coords(0.5, -0.1)

plt.tight_layout()


if len(argv)>1:
    filename = str(argv[1])
    plt.savefig(filename+'.pdf', format='pdf')

plt.show()
