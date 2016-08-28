from matplotlib import pyplot as plt
import argparse
import numpy as np
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-filename", type=str,
                    help="Name of file to read")
parser.add_argument("-save", action="store_true",
                    help="Save figure in figures folder")
parser.add_argument("-noPlot", action="store_true",
                    help="Do not show plot")
parser.add_argument("-start", type=float, default=0.5,
                    help="First x-value of the linear approximation is taken")
parser.add_argument("-stop", type=float, default=1,
                    help="Last x-value of the linear approximation is taken")
args = parser.parse_args()



infile = open(args.filename,"r")
infile.readline()
infile.readline()
content = infile.read().split()

timesteps = np.asarray([int(i)   for i in content[0::2]])
msd       = np.asarray([float(i) for i in content[1::2]])

timesteps = timesteps-timesteps[0]

dt   = 0.1    # Units though...
time = timesteps*dt
N    = len(time)
b    = int(np.floor(args.start * N))
e    = int(np.floor(args.stop  * N))

linearApproximation = np.polyfit(time[b:e], msd[b:e],1)
linearApproximation = np.poly1d(linearApproximation)

diffusionConstant = linearApproximation[1]/6.0



# Plot the mean square distance as a function of time.
if (parser.parse_args().save or not parser.parse_args().noPlot):
    figureLabel = "Linear approximation\n%s"%str(linearApproximation)
    plt.plot(time, msd, color="#36308A", linewidth=4)
    plt.xlabel('Time [fs]')
    plt.ylabel('Mean square displacement [AA]')
    plt.grid('on')
    plt.hold('on')
    plt.plot(time[b:e], np.polyval(linearApproximation, time[b:e]), color="red", linewidth=3, linestyle="--", label=figureLabel)
    plt.legend(loc=4)
    #plt.axis([dt, (timeSteps+1)*dt, 0, 1.1*max(msd)])
    plt.title("Diffusion constant: %.2g"%diffusionConstant)
    if parser.parse_args().save:
        filepath = "msdplotTEST"
        print "Figure saved as %s" %filepath
        plt.savefig(filepath)
    if not parser.parse_args().noPlot:
        plt.show()


print "Diffusion constant: ", diffusionConstant
