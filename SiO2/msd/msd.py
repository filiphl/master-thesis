# Script to plot mean square displacement

from matplotlib import pyplot as plt


def getDisplacement(T):
    filename = "msd_%d.txt"%T
    infile   = open(filename, "r")
    data     = infile.readlines()[-1].split()
    return int(data[0]), float(data[1])


temp = []
disp = []
j = 0
for i in xrange(500, 2000, 100):
    t,d = getDisplacement(i)
    temp.append(500+j*100)
    j+=1
    disp.append(d)

plt.plot(temp, disp)
plt.savefig("msdplot.png")
plt.show()
