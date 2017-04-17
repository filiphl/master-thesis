import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from os import listdir, system
import argparse
import numpy as np
import sys
import re





def D(filepath, PLOT=False, SAVE=False):
    temperature = int(re.findall(r'T(\d*)\.', filepath)[0])
    infile = open(filepath,'r')
    infile.readline()
    infile.readline()
    content = infile.read().split()

    timesteps = np.asarray([int(i)   for i in content[0::2]])
    msd       = np.asarray([float(i) for i in content[1::2]])

    timesteps = timesteps-timesteps[0]

    dt   = 0.002    # ps
    time = timesteps*dt
    N    = len(time)
    b    = int(np.floor(0.5 * N))
    e    = int(np.floor(1   * N))

    linearApproximation = np.polyfit(time[b:e], msd[b:e],1)
    linearApproximation = np.poly1d(linearApproximation)

    diffusionConstant = linearApproximation[1]/6.0



    # Plot the mean square distance as a function of time.
    if (PLOT or SAVE):
        figureLabel = "Linear approximation\n%s"%str(linearApproximation)
        plt.plot(time, msd, color="#36308A", linewidth=4)
        plt.xlabel('Time [ps]')
        plt.ylabel('Mean square displacement [AA]')
        plt.grid('on')
        plt.hold('on')
        plt.plot(time[b:e], np.polyval(linearApproximation, time[b:e]), color="red", linewidth=3, linestyle="--", label=figureLabel)
        plt.legend(loc=4)
        #plt.axis([dt, (timeSteps+1)*dt, 0, 1.1*max(msd)])
        plt.title("Temperature: %s\nDiffusion constant: %.2g"%(temperature,diffusionConstant))
        if SAVE:
            filepath = "msdplot_%05d.pdf"%temperature
            print "Figure saved as %s" %filepath
            plt.savefig(filepath)
            plt.clf()
        if PLOT:
            plt.show()

    return temperature, diffusionConstant


# -----------------------------------------------------------------------------#
folderPath = './dataFiles/'
directory  = listdir(folderPath)

temperature = []
diffusionConstant = []
msd = {}
for f in directory:
    t,d = D(folderPath+f, 0, 0)
    msd[t] = d


for T in sorted(msd.keys()):
    temperature.append(int(T))
    diffusionConstant.append(float(msd[T]))

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig, ax = plt.subplots()

ax.plot(temperature, diffusionConstant,
'-h',
color="#478684",
linewidth=4,
markersize=9,
markerfacecolor='#3b617c',
fillstyle='full')

plt.xlabel(r'Temperature [K]', fontsize=20)
plt.ylabel(r'Diffusion coefficient', fontsize=20)
plt.grid('on')
plt.ylim([min(diffusionConstant)*1.1,max(diffusionConstant)*1.1])
yticks = ax.get_yticks()
ax.set_yticks(yticks[2:])
plt.tick_params(axis='both', which='major', labelsize=20)


ax.xaxis.set_label_coords(0.5, -0.1)
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.tight_layout()


if True:
    plt.savefig('msd.pdf')
if True:
    plt.show()

system("mv msd*.pdf figures")
