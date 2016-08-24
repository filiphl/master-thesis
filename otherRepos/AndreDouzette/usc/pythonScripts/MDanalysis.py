from __future__ import division;
import numpy as np;
# import matplotlib.pyplot as plt;
import os;
import time;
import multiprocessing as mp;
try:
	import matplotlib.pyplot as plt;
except ImportError:
	pass
#end

def readEnergyFile(dirname):
	infile = open(dirname + "/energy.dat");
	out = [];
	#Kinetic energy
	out.append(float(infile.readline()));
	#Potential Energy
	out.append(float(infile.readline()));
	#Temperature
	out.append(float(infile.readline()));
	#<r^2>
	r2 = infile.readline().split(' ');
	for i in range(0, 3):
		if(r2[i] == "nan" or r2[i] == "-nan"):
			r2[i] = "0";
		#end
	#end
	out.append([float(r2[0]), float(r2[1]), float(r2[2])]);
	infile.close();
	return out;
#end

def readSimulationSize(runName):
	infile = open(runName + "/init/data.dat");
	infile.readline();
	L = [float(x) for x in infile.readline().split(' ')];
	infile.close();
	return np.array(L);
#end

def readParticlePositionsAndTypes(dirname):
	infile = open(dirname + "/combined.xyz");
	#Number of particles
	N = int(infile.readline());
	#discard comment
	infile.readline();
	r = np.zeros([N, 3]);
	t = np.zeros(N);
	for line in infile:
		data = line.split(' ');
		#id
		i = int(data[7]);
		#x, y, z positions
		x = float(data[1]);
		y = float(data[2]);
		z = float(data[3]);
		r[i, 0] = x;
		r[i, 1] = y;
		r[i, 2] = z;
		t[i] = int(data[0]);
	#end
	infile.close();
	return r, t;
#end

def readParticlePositions(dirname):
	infile = open(dirname + "/combined.xyz");
	#Number of particles
	N = int(infile.readline());
	#discard comment
	infile.readline();
	r = np.zeros([N, 3]);
	for line in infile:
		data = line.split(' ');
		#id
		i = int(data[7]);
		#x, y, z positions
		x = float(data[1]);
		y = float(data[2]);
		z = float(data[3]);
		r[i, 0] = x;
		r[i, 1] = y;
		r[i, 2] = z;
	#end
	infile.close();
	return r;
#end

def readParticleData(dirname):
	infile = open(dirname + "/combined.xyz");
	#Number of particles
	N = int(infile.readline());
	#discard comment
	infile.readline();
	r = np.zeros([N, 3]);
	v = np.zeros([N, 3]);
	F = np.zeros([N, 3]);
	t = np.zeros(N);
	for line in infile:
		data = line.split(' ');
		#id
		i = int(data[7]);
		#x, y, z positions
		r[i, 0] = float(data[1]);
		r[i, 1] = float(data[2]);
		r[i, 2] = float(data[3]);
		v[i, 0] = float(data[4]);
		v[i, 1] = float(data[5]);
		v[i, 2] = float(data[6]);
		F[i, 0] = float(data[8]);
		F[i, 1] = float(data[9]);
		F[i, 2] = float(data[10]);
		t[i] = int(data[0]);
	#end
	infile.close();
	return r, t, v, F;
#end

def readNumberOfParticles(runName):
	N = np.zeros(3);
	dirname = runName + "/" + "00000";
	infile = open(dirname + "/combined.xyz");
	infile.readline();
	infile.readline();	
	for line in infile:
		data = line.split(' ');
		N[int(data[0])] += 1;
	#end
	infile.close();
	for i in range(0, 3):
		if(N[i] == 0):
			N[i] = -1;
		#end
	#end
	return N;
#end

def getAnalysis(runName, last):
	N = last + 1;
	K = np.zeros(N);
	V = np.zeros(N);
	T = np.zeros(N);
	r2 = np.zeros([N, 3]);
	for i in range(0, last + 1):
		data = readEnergyFile(runName + "/" + ("%05d"%i));
		K[i] = data[0];
		V[i] = data[1];
		T[i] = data[2];
		r2[i] = data[3];
	#end
	return K, V, T, r2;
#end

def vmdVisualize():
	#Visualizing system using vmd
	os.system("vmd -e import.tcl");
#end

def calculateDisplacement(runName, last):
	print "Calculating <r^2>";
	#Calculates particle mean value of <r^2>
	#!Need to check if combined file exists
	L = readSimulationSize(runName);
	#Get first state
	r0, t = readParticlePositionsAndTypes(runName + "/00000");
	N = len(r0[:, 0]);
	#number of types
	Nt = int(np.max(t)) + 1;
	#number of particles of each type
	nt = np.zeros(Nt);
	for i in range(0, N):
		nt[t[i]] += 1;
	#end
	#Setting up position arrays
	r2 = np.zeros([last + 1, Nt]);
	r = np.copy(r0);
	rprev = r0;
	for i in range(1, last + 1):
		rnew = readParticlePositions(runName + ("/%05d"%i));
		for p in range(0, N):
			for d in range(0, 3):
				m = np.abs(rnew[p, d] - rprev[p, d]);
				#Check if moved across boundary
				if(m > np.abs(rnew[p, d] - rprev[p, d] - L[d])):
					r[p, d] = r[p, d] + rnew[p, d] - rprev[p, d] - L[d];
				elif(m > np.abs(rnew[p, d] - rprev[p, d] + L[d])):
					r[p, d] = r[p, d] + rnew[p, d] - rprev[p, d] + L[d];
				else:
					r[p, d] = r[p, d] + rnew[p, d] - rprev[p, d];
				#end
			#end
			r2[i, t[p]] += np.sum((r[p, :] - r0[p, :])**2);
		#end
		for j in range(0, Nt):
			r2[i, j] = r2[i, j]/nt[j];
		#end
		rprev = rnew;
	#end
	return r2;
#end

def radialDistributionThread(rank, states, bins, runName):
	Lvalues = readSimulationSize(runName);
	Lx = Lvalues[0];
	Ly = Lvalues[1];
	Lz = Lvalues[2];
	L = np.min(Lvalues);
	dr = 0.5*L/(bins - 1);
	dri = 1/dr;
	g = np.zeros([3, 3, bins]);
	
	position, t = readParticlePositionsAndTypes(runName + ("/%05d"%0));
	N = len(position);
	
	for s in states:
		if(rank == 0):
			print "%03d/%03d"%(s, states[-1]);
		#end
		position = readParticlePositions(runName + ("/%05d"%s));

		for i in range(0, N):
			r0 = position[i];
			rdiff = position[0:i] - r0;
			x = rdiff[:, 0];
			y = rdiff[:, 1];
			z = rdiff[:, 2];
			x = np.min([np.abs(x), np.abs(x - Lx), np.abs(x + Lx)], axis = 0);
			y = np.min([np.abs(y), np.abs(y - Ly), np.abs(y + Ly)], axis = 0);
			z = np.min([np.abs(z), np.abs(z - Lz), np.abs(z + Lz)], axis = 0);
			r = np.sqrt(x*x + y*y + z*z);
			n = (r*dri).astype(int);
			for j in range(0, i):
				if(n[j] < bins):
					g[t[i], t[j], n[j]] += 1;
				#end
			#end
		#end
	#end
	return g;
#end

def calculateRadialDistribution(runName, bins, last, states):
	if(states == "all"):
		states = np.arange(0, last + 1);
	elif(states == "last"):
		states = [last];
	#end
	states = np.array(states);
	
	Lvalues = readSimulationSize(runName);
	L = np.min(Lvalues);
	dr = 0.5*L/(bins - 1);
	g = np.zeros([3, 3, bins]);
	N = readNumberOfParticles(runName);
	n = np.arange(0, bins);
	V = 4*np.pi*dr*dr*dr*(n*n + n + 1/3);
	
	Ns = len(states);
	Ncpu = mp.cpu_count();
	pool = mp.Pool(processes = Ncpu);
	first = 0;
	process = [];
	for i in range(0, Ncpu):
		last = first + (Ns//Ncpu + (i < Ns%Ncpu));
		process.append(pool.apply_async(radialDistributionThread, args = (i, np.copy(states[first:last]), bins, runName)));
		first = last;
	#end
	pool.close();
	pool.join();
	for p in process:
		out = p.get();
		g += out;
	#end
	g[0, 1] += g[1, 0];
	g[0, 2] += g[2, 0];
	g[1, 2] += g[2, 1];
	g[1, 0] = g[1, 0]*0;
	g[2, 0] = g[2, 0]*0;
	g[2, 1] = g[2, 1]*0;
	for i in range(0, 3):
		g[i, i] = 2*g[i, i];
		for j in range(0, 3):
			g[i, j] = g[i, j]/(N[i]*V);
		#end
	#end
	return np.arange(0, bins)*dr, g/Ns;
#end

def calculateBonds(r, idSi, idO, idH):
	#No bonds
	if(len(idO) == 0):
		return np.array([]), np.array([]), np.array([]);
	elif(len(idSi) + len(idH) == 0):
		return np.array([]), np.array([]), np.array([]);
	#end
	#3BF cutoff
	rcSi = 2.75;
	rcH = 1.55;
	#Max number of bonds
	Nb = 8;
	#Number of particles
	N = len(r);
	#List of bonds
	bonds = np.zeros([N, Nb]) - 1;
	distances = np.zeros([N, Nb, 3]);
	bondNumber = np.zeros(N);
	
	rSi = r[idSi];
	rO = r[idO];
	rH = r[idH];
	
	r2cSi = rcSi*rcSi;
	r2cH = rcH*rcH;
	cells = np.zeros(3, dtype = int);
	
	####################
	#Finding Si-O bonds#
	####################
	if(len(idSi) != 0):
		cSi = np.array(rSi//rcSi, dtype = int);
		cO = np.array(rO//rcSi, dtype = int);
		for i in range(0, 3):
			m1 = np.max(cSi[:, i]);
			m2 = np.max(cO[:, i]);
			cells[i] = int(np.max([m1, m2])) + 1;
		#end
		L = rcSi*cells;
		# exit();
		#Creating cell structure
		cellSi = [];
		cellO = [];
		for i in range(0, cells[0]):
			cellSi.append([]);
			cellO.append([]);
			for j in range(0, cells[1]):
				cellSi[i].append([]);
				cellO[i].append([]);
				for k in range(0, cells[2]):
					cellSi[i][j].append([]);
					cellO[i][j].append([]);
				#end
			#end
		#end
		#Sorting particles into cells
		for i in range(0, len(idSi)):
			cellSi[cSi[i, 0]][cSi[i, 1]][cSi[i, 2]].append(idSi[i]);
		#end
		for i in range(0, len(idO)):
			cellO[cO[i, 0]][cO[i, 1]][cO[i, 2]].append(idO[i]);
		#end
		
		def findDistance2Si(pO, pSi):
			xO = r[pO, 0];
			yO = r[pO, 1];
			zO = r[pO, 2];
			xSi = r[pSi, 0];
			ySi = r[pSi, 1];
			zSi = r[pSi, 2];
			dx = xSi - xO;
			dx = np.min([np.abs(dx), np.abs(dx + L[0]), np.abs(dx - L[0])]);
			dy = ySi - yO;
			dy = np.min([np.abs(dy), np.abs(dy + L[1]), np.abs(dy - L[1])]);
			dz = zSi - zO;
			dz = np.min([np.abs(dz), np.abs(dz + L[2]), np.abs(dz - L[2])]);
			return dx*dx + dy*dy + dz*dz;
		#end
		
		#Finding bonds, looping through all neighbour cells
		for i in range(0, cells[0]):
			for j in range(0, cells[1]):
				for k in range(0, cells[2]):
					for di in range(-1, 2):
						i2 = (i + di + cells[0])%cells[0];
						for dj in range(-1, 2):
							j2 = (j + dj + cells[1])%cells[1];
							for dk in range(-1, 2):
								k2 = (k + dk + cells[2])%cells[2];
								for pO in cellO[i][j][k]:
									rOLocal = r[pO] - rcSi*np.array([i, j, k]);
									for pSi in cellSi[i2][j2][k2]:
										if(findDistance2Si(pO, pSi) <= r2cSi):
											rSiLocal = r[pSi] - rcSi*np.array([i2, j2, k2]);
											dr = rSiLocal - rOLocal + rcSi*np.array([di, dj, dk]);
											nSi = bondNumber[pSi];
											nO = bondNumber[pO];
											bonds[pO, nO] = pSi;
											bonds[pSi, nSi] = pO;
											distances[pO, nO, :] = dr;
											distances[pSi, nSi, :] = -dr;
											bondNumber[pO] += 1;
											bondNumber[pSi] += 1;
										#end
									#end
								#end
							#end
						#end
					#end
				#end
			#end
		#end
	#end
	###################
	#Finding O-H bonds#
	###################
	if(len(idH) != 0):
		cH = np.array(rH//rcH, dtype = int);
		cO = np.array(rO//rcH, dtype = int);
		for i in range(0, 3):
			m1 = np.max(cH[:, i]);
			m2 = np.max(cO[:, i]);
			cells[i] = int(np.max([m1, m2])) + 1;
		#end
		L = rcH*cells;
		#Creating cell structure
		cellH = [];
		cellO = [];
		for i in range(0, cells[0]):
			cellH.append([]);
			cellO.append([]);
			for j in range(0, cells[1]):
				cellH[i].append([]);
				cellO[i].append([]);
				for k in range(0, cells[2]):
					cellH[i][j].append([]);
					cellO[i][j].append([]);
				#end
			#end
		#end
		#Sorting particles into cells
		for i in range(0, len(idH)):
			cellH[cH[i, 0]][cH[i, 1]][cH[i, 2]].append(idH[i]);
		#end
		for i in range(0, len(idO)):
			cellO[cO[i, 0]][cO[i, 1]][cO[i, 2]].append(idO[i]);
		#end
		
		def findDistance2H(pO, pH):
			xO = r[pO, 0];
			yO = r[pO, 1];
			zO = r[pO, 2];
			xH = r[pH, 0];
			yH = r[pH, 1];
			zH = r[pH, 2];
			dx = xH - xO;
			dx = np.min([np.abs(dx), np.abs(dx + L[0]), np.abs(dx - L[0])]);
			dy = yH - yO;
			dy = np.min([np.abs(dy), np.abs(dy + L[1]), np.abs(dy - L[1])]);
			dz = zH - zO;
			dz = np.min([np.abs(dz), np.abs(dz + L[2]), np.abs(dz - L[2])]);
			return dx*dx + dy*dy + dz*dz;
		#end
		
		for i in range(0, cells[0]):
			for j in range(0, cells[1]):
				for k in range(0, cells[2]):
					for di in range(-1, 2):
						i2 = (i + di + cells[0])%cells[0];
						for dj in range(-1, 2):
							j2 = (j + dj + cells[1])%cells[1];
							for dk in range(-1, 2):
								k2 = (k + dk + cells[2])%cells[2];
								for pO in cellO[i][j][k]:
									rOLocal = r[pO] - rcH*np.array([i, j, k]);
									nO = bondNumber[pO];
									for pH in cellH[i2][j2][k2]:
										if(findDistance2H(pO, pH) <= r2cH):
											rHLocal = r[pH] - rcH*np.array([i2, j2, k2]);
											dr = rHLocal - rOLocal + rcH*np.array([di, dj, dk]);
											nO = bondNumber[pO];
											nH = bondNumber[pH];
											bonds[pO, nO] = pH;
											bonds[pH, nH] = pO;
											distances[pO, nO, :] = dr;
											distances[pH, nH, :] = -dr;
											bondNumber[pO] += 1;
											bondNumber[pH] += 1;
										#end
									#end
								#end
							#end
						#end
					#end
				#end
			#end
		#end
	#end
	return bonds, distances, bondNumber;
#end

def angularDistribution(bins, runName, last, states):
	if(states == "all"):
		states = np.arange(0, last + 1);
	elif(states == "last"):
		states = [last];
	#end
	angles = np.zeros([3, bins]);
	r, types = readParticlePositionsAndTypes(runName + "/%05d"%0);
	idSi = [];
	idO = [];
	idH = [];
	for i in range(0, len(types)):
		if(types[i] == 0):
			idSi.append(i);
		elif(types[i] == 1):
			idO.append(i);
		else:
			idH.append(i);
		#end
	#end
	
	L = readSimulationSize(runName);
	
	#Angle data
	osoAngles = np.zeros(bins);
	sosAngles = np.zeros(bins);
	sohAngles = np.zeros(bins);
	hohAngles = np.zeros(bins);
	angles = np.pi*np.arange(bins)/(bins - 1);
	
	da = (bins - 1)/np.pi;
	
	for s in states:
		r = readParticlePositions(runName + "/%05d"%s);
		bonds, distances, bondNumber = calculateBonds(r, idSi, idO, idH);
		#Bonds of O-Si-O
		for pSi in idSi:
			b = bonds[pSi];
			dr = distances[pSi];
			bn = int(bondNumber[pSi]);
			for i1 in range(bn):
				b1 = int(b[i1]);
				rv1 = dr[i1];
				r1 = np.linalg.norm(rv1);
				for i2 in range(i1):
					b2 = int(b[i2]);
					rv2 = dr[i2];
					r2 = np.linalg.norm(rv2);
					c = np.dot(rv1, rv2)/(r1*r2);
					#Avoiding floating point error for c = +/-(1 + eps)
					if(np.abs(c) > 1):
						t = np.pi;
					else:
						t = np.arccos(c);
					#end
					#Both bonds needs to be oxygen, as Si do not bond with anything else
					if(types[b1] == 1 and types[b2] == 1):
						osoAngles[int(da*t)] += 1;
					#end
				#end
			#end
		#end
		for pO in idO:
			b = bonds[pO];
			dr = distances[pO];
			bn = int(bondNumber[pO]);
			for i1 in range(bn):
				b1 = int(b[i1]);
				rv1 = dr[i1];
				r1 = np.linalg.norm(rv1);
				for i2 in range(i1):
					b2 = int(b[i2]);
					rv2 = dr[i2];
					r2 = np.linalg.norm(rv2);
					c = np.dot(rv1, rv2)/(r1*r2);
					#Avoiding floating point error for c = +/-(1 + eps)
					if(np.abs(c) > 1):
						t = np.pi;
					else:
						t = np.arccos(c);
					#end
					#Multiple particles bonds with oxygen
					if(types[b1] == 0 and types[b2] == 0):
						sosAngles[int(da*t)] += 1;
					elif(types[b1] == 0 and types[b2] == 2):
						sohAngles[int(da*t)] += 1;
					elif(types[b1] == 2 and types[b2] == 0):
						sohAngles[int(da*t)] += 1;
					elif(types[b1] == 2 and types[b2] == 2):
						hohAngles[int(da*t)] += 1;
					#end
				#end
			#end
		#end
	#end
	
	#Normalizing angle distributions
	n = np.sum(osoAngles);
	if(n != 0):
		osoAngles = osoAngles*da/n;
	#end
	n = np.sum(sosAngles);
	if(n != 0):
		sosAngles = sosAngles*da/n;
	#end
	n = np.sum(sohAngles);
	if(n != 0):
		sohAngles = sohAngles*da/n;
	#end
	n = np.sum(hohAngles);
	if(n != 0):
		hohAngles = hohAngles*da/n;
	#end
	angleDistribution = np.zeros([4, bins]);
	angleDistribution[0] = osoAngles;
	angleDistribution[1] = sosAngles;
	angleDistribution[2] = sohAngles;
	angleDistribution[3] = hohAngles;
	return angles, angleDistribution;
#end

def atomicDensity(dimension, bins, runName, last, states):
	if(states == "all"):
		states = np.arange(0, last + 1);
	elif(states == "last"):
		states = [last];
	#end
	data = np.zeros([3, bins]);
	s = readSimulationSize(runName);
	for s in states:
		r, types = readParticlePositionAndTypes(runName + "/%05d"%s);
		r = r[:, dimension];
		for j in range(0, len(types)):
			t = types[j];
			data[t, int(r[t]/s*bins)] += 1;
		#end
	#end
	return data;
#end

def silanolGroups(runName, last, states):
	if(states == "all"):
		states = np.arange(0, last + 1);
	elif(states == "last"):
		states = [last];
	#end
	
	r, types = readParticlePositionsAndTypes(runName + "/%05d"%0);
	idSi = [];
	idO = [];
	idH = [];
	for i in range(0, len(types)):
		if(types[i] == 0):
			idSi.append(i);
		elif(types[i] == 1):
			idO.append(i);
		else:
			idH.append(i);
		#end
	#end
	
	N = np.zeros(len(states));
	for i in range(len(states)):
		s = states[i];
		r = readParticlePositions(runName + "/%05d"%s);
		bonds, distances, bondNumber = calculateBonds(r, idSi, idO, idH);
		for pO in idO:
			si = False;
			h = False;
			for j in range(int(bondNumber[pO])):
				p = bonds[pO, j];
				if(types[p] == 0):
					si = True;
					# print str(p) + " is Si";
				#end
				if(types[p] == 2):
					h = True;
					# print str(p) + " is H";
				#end
			#end
			N[i] += si and h;
		#end
		print str(s) + ": " + str(int(N[i]));
	#end
	return N;
#end

def readForces(dirname):
	files = os.listdir(dirname);
	mfiles = [f for f in files if f[0:17] == "interatomicForces"];
	p = []; q = []; F = [];
	for f in mfiles:
		f = open(dirname + "/" + f, 'r');
		for line in f:
			split = line.split(' ');
			p.append(int(split[0]));
			q.append(int(split[1]));
			F.append([float(split[2]), float(split[3]), float(split[4])]);
		#end
	#end
	return np.array(p), np.array(q), np.array(F);
#end

def readMarkFile(runName, mark):
	name = runName + "/init/mark_%d.dat"%mark;
	if(not os.path.exists(name)):
		print "mark %d does not exist in the state %s"%(mark, runName);
		exit();
	#end
	f = open(name);
	return np.array([int(line.strip()) for line in f]);
#end

def forceDistribution(runName, bins, sender, states, direction, last):
	if(states == "all"):
		states = np.arange(0, last + 1);
	elif(states == "last"):
		states = [last];
	#end
	sender = readMarkFile(runName, sender);
	direction = np.array(direction);
	d = direction/np.sqrt(np.sum(direction*direction));
	
	L = readSimulationSize(runName);
	r0 = L/2;
	
	R = np.min([L[0], L[1]]);
	
	R = 0.25*np.sum(L - np.sum(d*L)*d);
	R2 = R*R;
	binsize = R*R/bins;
	A = binsize*np.pi;
	
	# data = np.zeros([3, bins]);
	data = np.zeros([3, bins, bins]);
	
	Lx = L[0];
	Ly = L[1];
	dx = Lx/(bins - 0);
	dy = Ly/(bins - 0);
	
	m = [28.0855, 15.9994, 1.00794];
	
	print " ";
	for s in states:
		print "\033[1A\033[K%d/%d"%(s, states[-1]);
		#get particle positions and forces
		# r = readParticlePositions(runName + "/%05d"%s);
		r, tp, v, _ = readParticleData(runName + "/%05d"%s);
		p, q, F = readForces(runName + "/%05d"%s);
		# print len(p);
		
		# for q in sender:
		# 	rq = r[q];
		# 	nx = int(rq[0]/dx);
		# 	ny = int(rq[1]/dy);
		# 	data[0, nx, ny] += rq[2];
		# #end
		# continue;
		for i in range(len(p)):
			
			#Pressure distribution from all particles onto selected
			rp = r[int(p[i])];
			rq = r[int(q[i])];
			nx = int(rp[0]/dx);
			ny = int(rp[1]/dy);
			data[0, nx, ny] += np.sum(F[i]*(rq - rp));
			data[1, nx, ny] += 0.5*m[int(tp[p[i]])]*np.sum(v[p[i]]**2);
			continue;
			
			if(not q[i] in sender):
				continue;
			#end
			rp = r[p[i]];
			rq = r[q[i]];
			nx = int(rp[0]/dx);
			ny = int(rp[1]/dy);
			
			drv = rq - rp;
			#Not exaclt drv, as cant be negative
			drv = np.array([np.min([np.abs(drv[i]), np.abs(drv[i] + L[i]), np.abs(drv[i] - L[i])]) for i in range(3)]);
			dr = np.sqrt(np.sum(drv*drv));
			#Need a factor r as F is stored as F/r
			Fv = F[i]*dr;
			data[:, nx, ny] += Fv;
			continue;
			data[0, nx, ny] += np.abs(dr*np.sum(Fv*drv));
			# rr = rq - r0 - d*np.sum((rq - r0)*d);
			# if(np.sum(rr*rr) > R2):
			# 	continue;
			# #end
			rr = np.array([rq[0] - r0[0], rq[1] - r0[1], 0]);
			# n = int(np.sum(rr*rr)/binsize);
			
			continue;
			#Radial force
			rr = rr/np.sqrt(np.sum(rr*rr));
			fz = np.sum(d*Fv);
			Fz = d*fz;
			Fv = Fv - Fz;
			fr = np.sum(rr*Fv);
			Ft = Fv - fr*rr;
			ft = np.sqrt(np.sum(Ft*Ft));
			data[0, n] += fz;
			data[1, n] += fr;
			data[2, n] += ft;
		#end
	#end
	# return R*np.sqrt(np.arange(bins)), data/A;
	x = np.arange(0, bins)*dx;
	y = np.arange(0, bins)*dy;
	return x, y, data/len(states);
#end