infile = open("cif2unit.xyz", "r")
outfile = open("output.xyz", "w")
for line in infile:
    word = line.split()
    for i in xrange(len(word)):
        try:
            word[i] = str(float(word[i])*7.12)
        except:
            continue
    outfile.write("    ".join(word)+"\n")
outfile.close()
