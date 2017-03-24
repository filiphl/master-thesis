from matplotlib import pyplot as plt
import numpy as np

infile = open('test.txt', 'r')

T = []
P = []
for line in infile:
    try:
        T.append(float(line.split()[1]))
        P.append(float(line.split()[-2]))
    except:
        continue
sP = np.std(P)
aP = np.mean(P)

plt.plot(T,P, '*')
plt.axis([0, 1800, -300, 300])
plt.show()
