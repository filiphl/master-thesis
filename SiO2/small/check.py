from math import sqrt

infile = open('forces.txt', 'r')

for i in xrange(3):
    infile.readline()

fs = []
for line in infile:
    col = line.split()
    if len(col) > 3:
        if int(col[3]) != 0:
            #print col[3]
            fs.append([float(col[4]), float(col[5]), float(col[6])])

sf = 0
for f in fs:
    a = 0
    for i in xrange(3):
        a += f[i]**2
    sf += sqrt(a)
print sf/len(fs)
