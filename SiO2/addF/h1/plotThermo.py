#Store the columns from a Thermo output in a dictionary, and plot this as a
#function of time steps. Also print their mean, though this is not interesting
#for all measurements.



import numpy as np
from matplotlib import pyplot as plt
from sys import argv


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("filename", type=str,
                    help="Name of slurm file.")
parser.add_argument("-mean", action="store_true",
                    help="Add line indicating mean value in plots.")
parser.add_argument("-single", type=int,
                    help="Plot single column.")
args = parser.parse_args()



c  = ['#D7DFC0','#AECFA2','#82BC92','#5FA38E','#49848B','#3F5F7F','#383D65','#2C1E3E']
cc = ['#AF97A9','#E5B3BB','#FCB7AE','#F7B490','#E3AA75','#C49859','#937E47','#5B5728']
#-------------------- Store data from file in dictionary -----------------------
infile = open(str(argv[1]), 'r')

slurmData = {}
key = []
for line in infile:
    col = line.split()
    if str(col[0]) == 'Step':
        nCols = len(col)
        for i in xrange(nCols):
            slurmData[str(col[i])] = []
            key.append(str(col[i]))
        break


for line in infile:
    col = line.split()
    try:
        for i in xrange(nCols):
            slurmData[key[i]].append(float(col[i]))
    except:
        continue


#-------------------------- Print their mean value -----------------------------

print '\n{0:<16} : {1:<16}'.format('Measurement', 'Mean value')
print '---------------------------------'
for i in xrange(1,nCols):
    print '{0:<16} : {1:<16}'.format(key[i], np.mean(slurmData[key[i]]))

#------------------------- Plot values vs time step ----------------------------


if args.single:
    plt.plot(slurmData['Step'], np.asarray(slurmData[key[args.single+1]]), linewidth=3, color=c[4])
    plt.ylabel(key[args.single+1], fontsize=14)
    plt.grid('on')
    print 'max: %.03f'%max(slurmData[key[args.single+1]])
else:
    fig, ax = plt.subplots(nCols-1, 1, sharex=True, figsize=(10,3*(nCols-1)))
    for i in xrange(0,nCols-1):
        ax[i].plot(slurmData['Step'], slurmData[key[i+1]], color=c[-i-1], linewidth=3)
        if args.mean:
            mean = np.mean(slurmData[key[i+1]])
            ax[i].axhline(y=mean, ls='--', color=cc[-i-1], linewidth=3)
        ax[i].set_ylabel(key[i+1], fontsize=14)

        #if key[i+1] in ['c_comTop[3]', 'c_comSphere[1]']:
        #     ax[i].set_ylim([mean-2, mean+2])

        #ax[i].set(aspect=40)
        ax[i].locator_params(axis='y', nbins=4)
        ax[i].grid('on')


plt.xlabel('Time step', fontsize=14)
plt.tight_layout()

plt.show()
