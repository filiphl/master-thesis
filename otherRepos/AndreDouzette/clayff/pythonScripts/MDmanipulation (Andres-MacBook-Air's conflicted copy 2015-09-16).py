from __future__ import division;
import numpy as np;
import os;

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

def combine(firstFile, secondFile, savefile, firstMarkFile, secondMarkFile, markSaveFile):
	data1 = readfile(firstFile);
	data2 = readfile(secondFile);
	N1 = len(data1);
	N2 = len(data2);
	#Updating ids of second file
	for i in range(0, N2):
		line = data2[i];
		split = line.split(' ');
		split[7] = str(N1 + int(split[7]));
		data2[i] = ' '.join(split) + "\n";
	#end
	#Adding all particles to outfile
	outfile = open(savefile, 'w');
	outfile.write(str(N1 + N2) + "\n<comment>\n");
	for line in data1:
		outfile.write(line.rstrip() + "\n");
	#end
	for line in data2:
		outfile.write(line.rstrip() + "\n");
	#end
	outfile.close();
	#Combining marked files if exists
	if(not os.path.exists(firstMarkFile) and not os.path.exists(secondMarkFile)):
		return;
	#end
	markfile = open(markSaveFile, 'w');
	if(os.path.exists(firstMarkFile)):
		mark1 = open(firstMarkFile, 'r');
		line = line;
		for line in mark1:
			markfile.write(line.strip() + "\n");
		#end
	#end
	if(os.path.exists(secondMarkFile)):
		mark2 = open(secondMarkFile, 'r');
		for line in mark2:
			line = line;
			split = line.split(' ');
			split[0] = str(N1 + int(split[0]));
			markfile.write(' '.join(split).strip() + "\n");
		#end
	#end
	markfile.close();
#end

def makeLens(filename, savefile, r0, h, d, L, direction, invert):
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
	mark = np.ones(len(t))*(1 - invert);
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
	#end
	saveParticleData(savefile, r, v, t, mark);
	return int(np.sum(1 - mark));
#end

def freeze(filename, markFilename, r1, r2, v, types):
	types = np.array(types);
	rm = np.zeros(3);
	rM = np.zeros(3);
	for d in range(0, 3):
		rm[d] = np.min([r1[d], r2[d]]);
		rM[d] = np.max([r1[d], r2[d]]);
	#end
	r, t = readParticlePositionsAndTypes(filename);
	markFile = open(markFilename, 'w');
	for i in range(0, len(t)):
		if(t[i] in types):
			if(np.all(r[i] >= rm) and np.all(r[i] <= rM)):
				markFile.write(str(i) + " " + str(v[0]) + " " + str(v[1]) + " " + str(v[2]) + "\n");
			#end
		#end
	#end
	markFile.close();
#end

def selectBox(filename, markFilename, r1, r2, types):
	pass
#end