from __future__ import division;
import numpy as np;
import os;
import shutil;
from math import floor;
import multiprocessing as mp;
import MDdata as mdd; mdd = mdd.usc();

def makeDirs(last, path, step = 1, reset = False):
	if(reset and os.path.exists(path)):
		shutil.rmtree(path);
	#end
	if(not os.path.exists(path)):
		os.makedirs(path);
	#end
	for i in range(0, last + step, step):
		dirname = path + ("/%05d"%i);
		if(not os.path.exists(dirname)):
			os.makedirs(dirname);
		#end
	#end
#end

def split(path, p, s, dest = None, source = None, remove = 0):
	if(dest == None):
		dest = path + "/00000";
	#end
	if(source == None):
		source = path + "/00000/combined.xyz";
	#end
	if(not os.path.exists(dest)):
		os.makedirs(dest);
	#end
	
	#reading file into 2D list of strings
	def readFile(filename, separator):
		infile = open(filename, "r");
		out = [];
		for line in infile:
			out.append(line.split(separator));
		#end
		infile.close();
		return out;
	#end

	#Fix input parameters
	#Read datafile
	datafilename = path + "/init/data.dat";
	data = readFile(datafilename, " ");
	px = p[0];
	py = p[1];
	pz = p[2];
	sx = s[0];
	sy = s[1];
	sz = s[2];
	
	sourcefilename = source;
	destfolder = dest;

	#Creates destfolder if it does not exist
	if not os.path.exists(destfolder):
		os.makedirs(destfolder);
	#end

	#number of processors
	procs = px*py*pz;
	#process file names
	procFilename = [];
	for i in range(0, procs):
		procFilename.append(destfolder + "/" + ("%03d"%i) + ".xyz");
	#end

	#Make arrays to vectorize calculations
	s = np.array([sx, sy, sz]);
	p = np.array([px, py, pz]);
	#spacial size of processors
	L = s/p;

	#open files
	sourcefile = open(sourcefilename);

	#number of particles in each processor
	procn = np.zeros(procs);
	#throw away number of particles and comments, 2 first lines;
	sourcefile.readline();
	sourcefile.readline();

	#finds number of particles in each processor
	proc = np.zeros(3);
	for line in sourcefile:
		tmp = np.array(line.split(" "));
		#coordinates
		r = tmp[1:4];
		for d in range(0, 3):
			fr = float(r[d]);
			proc[d] = int(fr/L[d]);
		#end
		rank = int(proc[0] + p[0]*(proc[1] + p[1]*proc[2]));
		procn[rank] += 1;
	#end

	#Writes number of particles and a comment to the process files
	for i in range(0, procs):
		tmpfile = open(procFilename[i], "w");
		tmpfile.write(str(procn[i]) + "\n<comment>\n");
		tmpfile.close();
	#end

	sourcefile.close();
	sourcefile = open(sourcefilename);
	#throw away number of particles and comments, 2 first lines
	sourcefile.readline();
	sourcefile.readline();
	tmpR = np.zeros(3);

	#separates each particle to corresponding processor
	for line in sourcefile:
		tmp = line.split(" ");
		#type
		t = str(int(tmp[0]) - remove);
		#coordinates
		r = tmp[1:4];
		for d in range(0, 3):
			fr = float(r[d]);
			proc[d] = int(fr/L[d]);
			tmpR[d] = fr - proc[d]*L[d];
		#end
		#Convert to local coordinates
		r = ["", "", ""];
		for d in range(0, 3):
			r[d] = str(tmpR[d]);
		#end
		rank = int(proc[0] + p[0]*(proc[1] + p[1]*proc[2]));
		#the rest remains the same
		rest = tmp[4:];
		tmpfile = open(procFilename[rank], "a");
		tmpfile.write(" ".join([t, " ".join(r), " ".join(rest)]));
		tmpfile.close();
	#end

	#close sourcefile
	sourcefile.close();
#end

def combineThread(path, states):
	if(len(states) == 0):
		return;
	#end
	combineFilename = "combined.xyz";
	#Read processor geometry and simulation vlume from data file
	datafile = open(path + "/init/data.dat", 'r');
	processors = [int(p) for p in datafile.readline().split(' ')];
	size = [float(s) for s in datafile.readline().split(' ')];
	datafile.close();
	px = processors[0];
	py = processors[1];
	pz = processors[2];
	#Calculate Lx, Ly, Lz for each processor
	procs = px*py*pz;
	s = np.array(size);
	p = np.array(processors);
	#spacial size of processors
	L = s/p;
	#Coordinates of processors
	r0 = np.zeros([procs, 3]);
	for i in range(px):
		for j in range(py):
			for k in range(pz):
				rank = i + px*j + px*py*k;
				r0[rank, :] = np.array([i, j, k])*L;
			#end
		#end
	#end
	# for rank in range(0, procs):
	# 	#Finding x, y, z latice positions of processes
	# 	i = rank%p[0];
	# 	j = (rank//p[0])%p[1];
	# 	k = rank//(p[0]*p[1]);
	# 	lr = np.array([i, j, k]);
	# 	r0[rank, :] = lr*L;
	# #end
	#get number of particles
	datafile = open(path + "/init/data.dat");
	lines = [l for l in datafile];
	N = int(lines[4]);
	datafile.close();
	particles = ['']*N;
	for s in states:
		stateDir = path + "/%05d"%s;
		combineFile = stateDir + "/" + combineFilename;
		if(os.path.exists(combineFile)):
			continue;
		#end
		#Read particle data from processor files
		for rank in range(0, procs):
			particleFile = open(stateDir + "/%03d.xyz"%rank, 'r');
			#Discard two first lines
			particleFile.readline();
			particleFile.readline();
			for line in particleFile:
				line = line.strip();
				if(line == ""):
					continue;
				#end
				split = line.split(' ');
				i = int(split[7]);
				#Transform processor coordinates to global coordinates
				split[1] = str(r0[rank, 0] + float(split[1]));
				split[2] = str(r0[rank, 1] + float(split[2]));
				split[3] = str(r0[rank, 2] + float(split[3]));
				particles[i] = ' '.join(split);
			#end
			particleFile.close();
		#end
		
		#Write particle data to a combined file
		combineFile = open(combineFile, 'w');
		combineFile.write(str(N) + "\n");
		combineFile.write("<comment>\n");
		for i in range(0, N):
			combineFile.write(particles[i].strip() + "\n");
		#end
		combineFile.close();
	#end
#end

#Update to this!
def combine(path, last, states = "all", processors = "all"):
	if(states == "all"):
		states = range(0, last + 1);
	elif(states == "last"):
		states = [last];
	#end
	if(processors == "all"):
		Ncpu = mp.cpu_count();
	else:
		Ncpu = int(processors);
	#end
	
	#Combines files using a parallell procedure
	Ns = len(states);
	jobs = []
	first = 0;
	for i in range(Ncpu):
		last = first + (Ns//Ncpu + (i < Ns%Ncpu));
		p = mp.Process(target = combineThread, args = (path, np.copy(states[first:last])));
		jobs.append(p)
		p.start()
		first = last;
	#end
	for j in jobs:
		j.join();
	#end
#end

def makeVmdData(path, last, step = 1):
	infile = "combined.xyz";
	filename = "./vmd.dat";
	dfile = open(filename, "w");
	dfile.write(path + "\n");
	dfile.write(infile + "\n");
	dfile.write(str(last) + "\n");
	dfile.write(str(step));
#end




