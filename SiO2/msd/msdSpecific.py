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
        if linearApproximation[0] < 0:
            figureLabel = "Linear approximation\n%.5ft - %.5f"%(linearApproximation[1], abs(linearApproximation[0]))
        else:
            figureLabel = "Linear approximation\n%.5ft + %.5f"%(linearApproximation[1], abs(linearApproximation[0]))
        clr  = ['#D7DFC0','#AECFA2','#82BC92','#5FA38E','#49848B','#3F5F7F','#383D65','#2C1E3E']
        fig, ax = plt.subplots()
        plt.plot(time, msd, color=clr[5], linewidth=4)
        plt.xlabel(r'Time [ps]', fontsize=20)
        plt.ylabel(r'Mean square displacement [\AA]', fontsize=20)
        plt.grid(1)
        plt.hold(1)
        plt.plot(time[b:e], np.polyval(linearApproximation, time[b:e]), color=[0.7, 0.2, 0.2], linewidth=3, linestyle="--")#, label=figureLabel)
        #plt.legend(loc=2, fontsize=15)
        plt.tick_params(axis='both', which='major', labelsize=20)
        plt.ylim([0, 5])
        plt.xlim([0, 40])
        #plt.axis([dt, (timeSteps+1)*dt, 0, 1.1*max(msd)])
        plt.title("{0:<20s}{1:<10d}K\n{2:<20s}{3:<10.3f}".format('Temperature:    ', temperature, 'Diffusion constant:', diffusionConstant), fontsize=20)
        ax.xaxis.set_label_coords(0.5, -0.13)
        ax.yaxis.set_label_coords(-0.08, 0.5)
        plt.tight_layout()
        if SAVE:
            filepath = "msdplot_%05d.pdf"%temperature
            print "Figure saved as %s" %filepath
            plt.savefig(filepath, format='pdf', transparent=True)
            plt.close()
        if PLOT:
            plt.show()
            plt.close()

    return temperature, diffusionConstant


# -----------------------------------------------------------------------------#
clr = ['#D7DFC0','#AECFA2','#82BC92','#5FA38E','#49848B','#3F5F7F','#383D65','#2C1E3E']

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

folderPath = './dataFiles/'
directory  = listdir(folderPath)

temperature = []
diffusionConstant = []
msd = {}
c = 0
for f in directory:
    if c>-1:
        t,d = D(folderPath+f, 0, 1)
        msd[t] = d
    c+=1

for T in sorted(msd.keys()):
    temperature.append(int(T))
    diffusionConstant.append(float(msd[T]))

fig, ax = plt.subplots()

ax.plot(temperature, diffusionConstant,
'-h',
color=clr[5],
linewidth=4,
markersize=9,
markerfacecolor=clr[6],
fillstyle='full')

plt.xlabel(r'Temperature [K]', fontsize=20)
plt.ylabel(r'Diffusion coefficient', fontsize=20)
plt.grid('on')
plt.ylim([min(diffusionConstant)*1.1,max(diffusionConstant)*1.1])
plt.xlim([1500, 2000])
yticks = ax.get_yticks()
ax.set_yticks(yticks[2:])
plt.tick_params(axis='both', which='major', labelsize=20)


ax.xaxis.set_label_coords(0.5, -0.13)
ax.yaxis.set_label_coords(-0.15, 0.5)
plt.tight_layout()


if True:
    plt.savefig('msdTransparent.pdf', format='pdf', transparent=True)
if True:
    plt.show()

system("mv msd*.pdf figures/")
