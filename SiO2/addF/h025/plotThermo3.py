#Store the columns from a Thermo output in a dictionary, and plot this as a
#function of time steps. Also print their mean, though this is not interesting
#for all measurements.



import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from sys import argv
import re
import os

'''

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("filename", type=str,
                    help="Name of slurm file.")
parser.add_argument("-mean", action="store_true",
                    help="Add line indicating mean value in plots.")
parser.add_argument("-single", type=int,
                    help="Plot single column.")
args = parser.parse_args()
'''


c  = ['#D7DFC0','#AECFA2','#82BC92','#5FA38E','#49848B','#3F5F7F','#383D65','#2C1E3E']
cc = ['#AF97A9','#E5B3BB','#FCB7AE','#F7B490','#E3AA75','#C49859','#937E47','#5B5728']
#-------------------- Store data from file in dictionary -----------------------

def load(filepath):

    slurmData = {}
    key = []

    infile = open(filepath, 'r')
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

    infile.close()
    content = open(filepath, 'r').read()
    v = float(re.findall(r'vel=-([^\s]*)', content)[0])

    return slurmData, key, v

#-------------------------- Print their mean value -----------------------------
def printData(slurmData, key, v):
    Fmax = max(slurmData['f_smdID[1]'])
    N = np.mean(slurmData[key[-1]])
    u=Fmax/N
    print '\n{0:<16}{1:<16}{2:<16}{3:<16}{4:<16}'.format('Step', 'velocity', 'Fmax', 'N', 'u')
    print   '{0:<16.0f}{1:<16.1f}{2:<16.2f}{3:<16.2f}{4:<16.2f}'.format(slurmData['Step'][0], v, Fmax, N, u)

#------------------------- Plot values vs time step ----------------------------


def scaleF(F, v):
    F /= 0.0702*v + 0.3992
    return F





def plotV(slurmData, v, ci):
    t = (np.asarray(slurmData['Step'])-float(min(slurmData['Step'])))*0.002 #[ps]
    d = v*t #AA
    F = np.asarray(slurmData['f_smdID[1]'])
    N = np.mean(slurmData[key[-1]])
    plt.plot(scaleF(d,v), scaleF(F,v)/N, linewidth=3, label='$v=%d$'%v, color=c[ci])
    plt.ylabel('$\\frac{F}{N\\mu(v)}$', rotation=0, fontsize=28)#'Force [eV/\AA]', fontsize=16)
    plt.xlabel('$\\frac{t v}{\\mu(v)}$', fontsize=28)
    plt.grid('on')
    plt.tick_params(axis='both', which='major', labelsize=16)

    ax = plt.gca()
    ax.yaxis.set_label_coords(-0.12, 0.455)



if __name__ == '__main__':
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    dirpath = 'sameNdifferentV/'
    directory = os.listdir('sameNdifferentV/')
    ci = 2
    for filepath in sorted(directory):
        slurmData, key, v = load(dirpath+filepath)
        printData(slurmData, key, v)
        plotV(slurmData, v, ci)
        plt.hold('on')
        ci+=1
    plt.legend()
    plt.tight_layout()
    plt.show()
