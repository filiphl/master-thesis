import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc

def LJ(r, theta, epsilon):
  return 4*epsilon*( (theta/r)**12 - (theta/r)**6 )

r = np.linspace(0.95,3,1000)
V = LJ(r,1,1)

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig, ax = plt.subplots(figsize=(15,5))
ax.plot(r,V,linewidth=5, color='#49848B')
ax.grid('off')
plt.xlabel(r'$r$', fontsize=35)
plt.ylabel(r'$V(r)$', fontsize=30)
ax.set_xticks(range(4))
ax.set_yticks([0])
ax.set_yticklabels([0], fontsize=25)
ax.set_xticklabels(range(4), fontsize=25)
plt.axis([0,3,-1.15, 1.15])
fig.tight_layout()
plt.savefig('LJpython.pdf', transparency=True)
plt.show()
