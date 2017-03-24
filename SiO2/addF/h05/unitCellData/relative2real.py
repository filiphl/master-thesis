import sys

infile = open('unit.xyz', 'r')
L = float(sys.argv[1])
outfile = open('unit%s.xyz'%str(L), 'w')
for line in infile:
    col = line.split()
    try:
        out = '$atom:{0:<10}@atom:{1:<10}{2:<10}{3:<10}{4:<10}\n'.format(col[0], col[0], float(col[1])*L, float(col[2])*L, float(col[3])*L)
    except:
        out = line
    outfile.write(out)
