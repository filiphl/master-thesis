from __future__ import division;
import numpy as np;
import os;
import shutil;

#function to calculate bonds
#function to make a square hole
#function to make a cylindrical hole
#function to make a spherical hole
#function to make inverted holes
#functions to fill in OH groups at SiO2 surface
#function to remove water molecule if parts has already been removed

def readParticlePositionsAndTypes(filename):
	infile = open(filename);
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
		#Remove 1, as 1 is added due to vmd
		t[i] = int(data[0]);
	#end
	infile.close();
	return r, t;
#end

def readParticlePositions(filename):
	infile = open(filename);
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

def saveParticleData(savefile, r, v, t, mark):
	id = 0;
	sfile = open(savefile, 'w');
	sfile.write(str(int(np.sum(1 - mark))) + "\n");
	sfile.write("<Comment>\n");
	for i in range(0, len(t)):
		if(not mark[i]):
			x = str(r[i, 0]);
			y = str(r[i, 1]);
			z = str(r[i, 2]);
			vx = str(v[i, 0]);
			vy = str(v[i, 1]);
			vz = str(v[i, 2]);
			sfile.write(str(int(t[i])) + " " + x + " " + y + " " + z + " " + vx + " " + vy + " " + vz + " " + str(id) + "\n");
			id += 1;
		#end
	#end
	sfile.close();
#end

def readParticleData(filename):
	infile = open(filename);
	#Number of particles
	N = int(infile.readline());
	#discard comment
	infile.readline();
	r = np.zeros([N, 3]);
	v = np.zeros([N, 3]);
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
		t[i] = int(data[0]);
	#end
	infile.close();
	return r, v, t;
#end

def readfile(filename):
	infile = open(filename, "r");
	#Discard two first lines
	infile.readline();
	infile.readline();
	#Adding all other lines to list
	out = [];
	for line in infile:
		out.append(line);
	#end
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

def moveParticles(loadfile, savefile, dr, L):
	infile = open(loadfile, 'r');
	outfile = open(savefile, 'w');
	outfile.write(infile.readline());
	outfile.write(infile.readline());
	for line in infile:
		split = line.split(' ');
		x = float(split[1]);
		y = float(split[2]);
		z = float(split[3]);
		x = str((x + L[0] + dr[0])%L[0]);
		y = str((y + L[1] + dr[1])%L[1]);
		z = str((z + L[2] + dr[2])%L[2]);
		out = split[0] + " " + x + " " + y + " " + z;
		for i in range(4, len(split)):
			out = out + " " + split[i];
		#end
		outfile.write(out);
	#end
	infile.close();
	outfile.close();
#end

def combine(load1, load2, save):
	lastState1 = len(os.listdir(load1)) - 2;
	lastState2 = len(os.listdir(load2)) - 2;
	data1 = readfile(load1 + "/%05d/combined.xyz"%lastState1);
	data2 = readfile(load2 + "/%05d/combined.xyz"%lastState1);
	N1 = len(data1);
	N2 = len(data2);
	if(not os.path.exists(save + "/00000")):
		os.makedirs(save + "/00000");
	#end
	if(not os.path.exists(save + "/init")):
		os.makedirs(save + "/init");
	#end
	#Updating ids of second file
	for i in range(0, N2):
		line = data2[i];
		split = line.split(' ');
		split[7] = str(N1 + int(split[7]));
		data2[i] = ' '.join(split);
	#end
	#Adding all particles to outfile
	outfile = open(save + "/00000/combined.xyz", 'w');
	outfile.write(str(N1 + N2) + "\n<comment>\n");
	for line in data1:
		outfile.write(line.strip() + "\n");
	#end
	for line in data2:
		outfile.write(line.strip() + "\n");
	#end
	outfile.close();
	#Combining marked files if exists
	initfiles1 = os.listdir(load1 + "/init");
	initfiles2 = os.listdir(load2 + "/init");
	markfiles1 = [];
	markfiles2 = [];
	for i in range(len(initfiles1)):
		if(initfiles1[i][0:4] == "mark"):
			markfiles1.append(initfiles1[i]);
		#end
	#end
	for i in range(len(initfiles2)):
		if(initfiles2[i][0:4] == "mark"):
			markfiles2.append(initfiles2[i]);
		#end
	#end
	#Copy mark files of first state
	for mark in markfiles1:
		shutil.copyfile(load1 + "/init/" + mark, save + "/init/" + mark);
	#end
	#Adjust indexes of markfiles of second state
	for mark in markfiles2:
		fin = open(load2 + "/init/" + mark, 'r');
		fout = open(save + "/init/" + mark, 'w');
		for line in fin:
			fout.write(str(int(line.strip()) + N1) + "\n");
		#end
		fin.close();
		fout.close();
	#end
#end

def makeLens(load, save, r0, h, d, L, direction, invert):
	if(not os.path.exists(save + "/00000")):
		os.makedirs(save + "/00000");
	#end
	if(not os.path.exists(save + "/init")):
		os.makedirs(save + "/init");
	#end
	lastState = len(os.listdir(load)) - 2;
	filename = load + "/%05d/combined.xyz"%lastState;
	savefile = save + "/00000/combined.xyz";
	r, v, t = readParticleData(filename);
	r0 = np.array(r0);
	direction = np.array(direction);
	direction = direction/np.linalg.norm(direction);
	#Distance between two circles
	D = (h*h - d*d)/(4*d);
	#Radius of circles
	R = (h*h + d*d)/(4*d);
	r1 = r0 - D*direction;
	r2 = r0 + D*direction;
	R2 = R*R;
	#Marked for deletion
	mark = np.ones(len(t))*(1 - invert);
	id = 0;
	newp = np.zeros(len(t));
	for i in range(0, len(t)):
		dx1 = r[i, 0] - r1[0];
		dy1 = r[i, 1] - r1[1];
		dz1 = r[i, 2] - r1[2];
		dx2 = r[i, 0] - r2[0];
		dy2 = r[i, 1] - r2[1];
		dz2 = r[i, 2] - r2[2];
		if(dx1*dx1 + dy1*dy1 + dz1*dz1 < R2 and dx2*dx2 + dy2*dy2 + dz2*dz2 < R2):
			mark[i] = invert;
		#end
		if(not mark[i]):
			newp[i] = id;
			id += 1;
		else:
			newp[i] = -1;
		#end
	#end
	saveParticleData(savefile, r, v, t, mark);
	#Adjust indexes of markfiles
	initfiles = os.listdir(load + "/init");
	markfiles = [];
	for i in range(len(initfiles)):
		if(initfiles[i][0:4] == "mark"):
			markfiles.append(initfiles[i]);
		#end
	#end
	for m in markfiles:
		fin = open(load + "/init/" + m, 'r');
		fout = open(save + "/init/" + m, 'w');
		for line in fin:
			p = int(line.strip());
			#If not marked for deletion
			if(newp[p] != -1):
				fout.write(str(int(newp[p])) + "\n");
			#end
		#end
		fin.close();
		fout.close();
	#end
	return int(np.sum(1 - mark));
#end

def displaceRadially(load, save, R, direction):
	#Normalizing direction
	direction = np.array(direction);
	d = direction/np.sqrt(np.sum(direction*direction));
	R2 = R*R;
	lastState = len(os.listdir(load)) - 2;
	loadfile = load + "/%05d/combined.xyz"%lastState;
	savefile = save + "/00000/combined.xyz";
	infile = open(loadfile, 'r');
	outfile = open(savefile, 'w');
	outfile.write(infile.readline());
	outfile.write(infile.readline());
	
	L = readSimulationSize(load);
	r0 = L/2;
	for line in infile:
		split = line.split(' ');
		x = float(split[1]);
		y = float(split[2]);
		z = float(split[3]);
		#Displacement 
		r = np.array([x, y, z]);
		a = r - r0;
		r2 = np.sum(a*a) - np.sum(a*d)**2;
		if(r2 > R2):
			h = 0;
		else:
			h = R - np.sqrt(R2 - r2);
		#end
		#Getting displacement vector
		dr = h*d;
		x = str((x + L[0] + dr[0])%L[0]);
		y = str((y + L[1] + dr[1])%L[1]);
		z = str((z + L[2] + dr[2])%L[2]);
		out = split[0] + " " + x + " " + y + " " + z;
		for i in range(4, len(split)):
			out = out + " " + split[i];
		#end
		outfile.write(out);
	#end
	infile.close();
	outfile.close();
#end 

def selectBox(filename, markFilename, r1, r2, types):
	if(types == "all"):
		types = [0, 1, 2];
	#end
	rm = np.array([np.min([r1[d], r2[d]]) for d in range(3)])
	rM = np.array([np.max([r1[d], r2[d]]) for d in range(3)]);
	r, t, = readParticlePositionsAndTypes(filename);
	markFile = open(markFilename, 'w');
	for i in range(0, len(t)):
		if(t[i] in types):
			if(np.all(r[i] >= rm) and np.all(r[i] <= rM)):
				markFile.write(str(i) + "\n");
			#end
		#end
	#end
	markFile.close();
#end

def freeze(markFilename, frozenFilename, v):
	#Open existing file if there is any, else create new.
	if(os.path.exists(frozenFilename)):
		frozenFile = open(frozenFilename, 'a');
	else:
		frozenFile = open(frozenFilename, 'w');
	#end
	markFile = open(markFilename, 'r');
	for line in markFile:
		p = int(line.strip());
		frozenFile.write("%d %f %f %f\n"%(p, v[0], v[1], v[2]));
	#end
	markFile.close();
	frozenFile.close();
#end