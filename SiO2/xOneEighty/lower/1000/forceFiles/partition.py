from sys import argv

filename = str(argv[1])

infile = open(filename, 'r')

header = ''
for line in range(3):
    header += infile.readline()	# Dump file


for line in infile:
    col = line.split()
    if len(col)<4:
        try:
            timestep = int(col[0])
            outname = 'forces'+str(timestep)+'.txt'
            outfile = open(outname, 'w')
            outfile.write(header)
        except:
            print 'RIGHT HERE!', line
            continue
    outfile.write(line)

