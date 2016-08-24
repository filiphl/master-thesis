from __future__ import division;
import numpy as np;
import os;
import random;
from math import pi as tau; tau = 2*tau;
from math import cos, sin;

#!Include removal of total velocity

##private
def createDir(dirname):
	if(not os.path.exists(dirname)):
		os.makedirs(dirname);
	#end
#end

#Removes the velocity of mass the center
def removeVelocity(filename):
	infile = open(filename);
	N = int(infile.readline());
	comment = infile.readline();
	v = np.zeros([N, 3]);
	lines = [];
	t = np.zeros(N);
	p = 0;
	#Finding velocities
	for line in infile:
		lines.append(line);
		data = line.split(' ');
		t[p] = int(data[0]);
		v[p, 0] = float(data[4]);
		v[p, 1] = float(data[5]);
		v[p, 2] = float(data[6]);
		p += 1;
	#end
	infile.close();
	#Finding mean values
	Nt = int(np.max(t)) + 1;
	vx = np.zeros(Nt);
	vy = np.zeros(Nt);
	vz = np.zeros(Nt);
	nt = np.zeros(Nt);
	for p in range(0, N):
		vx[t[p]] += v[p, 0];
		vy[t[p]] += v[p, 1];
		vz[t[p]] += v[p, 2];
		nt[t[p]] += 1;
	#end
	vx = vx/nt;
	vy = vy/nt;
	vz = vz/nt;
	#Removing mean velocity
	outfile = open(filename, 'w');
	outfile.write(str(N) + "\n");
	outfile.write(comment);
	for p in range(0, N):
		v[p, 0] = v[p, 0] - vx[t[p]];
		v[p, 1] = v[p, 1] - vy[t[p]];
		v[p, 2] = v[p, 2] - vz[t[p]];
		data = lines[p].split(' ');
		line = data[0] + ' ';
		line = line + data[1] + ' ';
		line = line + data[2] + ' ';
		line = line + data[3] + ' ';
		line = line + str(v[p, 0]) + ' ';
		line = line + str(v[p, 1]) + ' ';
		line = line + str(v[p, 2]);
		for i in range(7, len(data)):
			line = line + ' ' + data[i];
		#end
		outfile.write(line);
	#end
	outfile.close();
#end

def water(s, c, path, T, border):
	typeO = 3;
	typeH = 2;
	if(len(s) != 3 or len(c) != 3):
		print "Error! s, and c must be 3d vectors";
		return;
	#end
	s = np.array(s);
	c = np.array(c);
	#size of simulation space
	sx = s[0];
	sy = s[1];
	sz = s[2];
	#number of atoms in each direction
	cx = int(c[0]);
	cy = int(c[1]);
	cz = int(c[2]);
	#size of cells
	L = s/c;
	#If packed too tight together
	for l in L:
		if(l < 1.633):
			print L;
			print "Molecules will overlap. Maximum density is 6.872g/cm.";
			exit();
		#end
	#end
	#rcut in terms of R(LJ)
	rcut = 3;
	#Number of particle types
	numTypes = 2;
	
	MO = 15.999;
	MH = 1.008;
	
	#std of velocity distribution
	sigma = 0.000911837*np.sqrt(T/(MO + 2*MH));
	
	N = 3*cx*cy*cz;

	#-----------#
	# Functions #
	#-----------#
	def waterMolecule(B, border):
		def rotate(r, t):
			tx = 2*np.pi*t[0];
			ty = np.arcsin(2*t[1] - 1);
			tz = 2*np.pi*t[2];
			x = r[0]; y = r[1]; z = r[2];
			#Rotating about x axis
			c = np.cos(tx); s = np.sin(tx);
			tmp = y*c - z*s;
			z   = z*c + y*s;
			y   = tmp;
			#Rotating about y axis
			c = np.cos(ty); s = np.sin(ty);
			tmp = x*c + z*s;
			z   = z*c - x*s;
			x   = tmp;
			#Rotating about z axis
			c = np.cos(tz); s = np.sin(tz);
			tmp = x*c - y*s;
			y   = y*c + x*s;
			x   = tmp;
			return np.array([x, y, z]);
		#end
		pi = np.pi;
		d = 0.9584; #Distance between O and H
		theta = 104.45/180*np.pi; #H-O-H angle
		c = np.cos(0.5*theta);
		s = np.sin(0.5*theta);
		x0 = 0.5*d/c;
		y0 = 0;
		x1 = x0 - d*c;
		x2 = x1;
		y1 = y0 + d*s;
		y2 = y0 - d*s;
		r0 = np.array([x0, y0, 0]);
		r1 = np.array([x1, y1, 0]);
		r2 = np.array([x2, y2, 0]);
		rotation = np.random.random(size = 3);
		r0 = rotate(r0, rotation);
		r1 = rotate(r1, rotation);
		r2 = rotate(r2, rotation);
		rm = np.array([np.min([r0[d], r1[d], r2[d]]) for d in range(3)]);
		rM = np.array([np.max([r0[d], r1[d], r2[d]]) for d in range(3)]);
		r0 = r0 - rm;
		r1 = r1 - rm;
		r2 = r2 - rm;
		L = rM - rm;
		b = 0.5*(B - L)*border;
		displacement = b + (B - L - 2*b)*np.random.random(size = 3);
		r0 = r0 + displacement;
		r1 = r1 + displacement;
		r2 = r2 + displacement;
		return np.array([r0, r1, r2]);
	#end
	
	def randomV(sigma):
		return random.gauss(0, sigma);
	#end
	
	def randomVelocity(sigma):
		return np.array([randomV(sigma), randomV(sigma), randomV(sigma)]);
	#end

	def saveToFile(r, v, ID, type, bonds, ifile):
		ifile.write(str(type) + " ");
		for i in range(0, 3):
			ifile.write(str(r[i]) + " ");
		#end
		for i in range(0, 3):
			ifile.write(str(v[i]) + " ");
		#end
		ifile.write(str(ID));
		#Writing bonds
		for i in range(0, len(bonds)):
			ifile.write(" " + str(bonds[i]));
		#end
		ifile.write("\n");
	#end

	#------#
	# main #
	#------#
	#Creating folder and files
	createDir(path);
	foldername = path + "/init";
	createDir(foldername);
	createDir(path + "/00000");
	initfilename = path + "/00000/combined.xyz";
	
	bondFile = open(foldername + "/bonds.dat", 'w');
	
	initfile = open(initfilename, 'w');
	initfile.write(str(N) + "\ntype x y z vx vy vz ID {bonds}\n");
	
	for i in range(0, 3):
		if(L[i] < 0):
			print("Cell sizes are too small.");
			#Quits program
			exit();
		#end
	#end
	
	#create file and directory
	for i in range(0, cx):
		for j in range(0, cy):
			for k in range(0, cz):
				cell = i + cx*(j + cy*k);
				R = L*[i, j, k];
				
				r0, r1, r2 = waterMolecule(L, border);
				
				r0 = R + r0;
				r1 = R + r1;
				r2 = R + r2;
				
				cell = 3*cell;
				v = randomVelocity(sigma);
				saveToFile(r0, v, cell, typeO, [cell + 1, cell + 2], initfile);
				saveToFile(r1, v, cell + 1, typeH, [cell], initfile);
				saveToFile(r2, v, cell + 2, typeH, [cell], initfile);
				bondFile.write(str(cell) + " " + str(cell + 1) + " " + str(cell + 2) + "\n");
				bondFile.write(str(cell + 1) + " " + str(cell) + "\n");
				bondFile.write(str(cell + 2) + " " + str(cell) + "\n");
			#end
		#end
	#end
	initfile.close();
	bondFile.close();
	return N;
#end

def cristobalite(s, c, T, path):
	typeSi = 0;
	typeO = 1;
	if(len(s) != 3 or len(c) != 3):
		return;
	#end
	#------------#
	# input data #
	#------------#
	#number of cells
	cx = int(c[0]);
	cy = int(c[1]);
	cz = int(c[2]);

	#-------------------#
	# derived constants #
	#-------------------#
	#size of simulation space
	b = s/c;
	N = 24*cx*cy*cz;
	#Import from data file/class?
	MSi = 28.0855;
	MO = 15.9994;
	sigmaSi = 0.000911837*np.sqrt(T/MSi);
	sigmaO = 0.000911837*np.sqrt(T/MO);
	
	#-----------#
	# functions #
	#-----------#
	#!This is for gasses, not crystals, change this?
	def randomV(sigma):
		return random.gauss(0, sigma);
	#end
	
	def randomVelocity(sigma):
		return np.array([randomV(sigma), randomV(sigma), randomV(sigma)]);
	#end
	
	def unitCell(i, j, k):
		#Cell indexes to the right
		pi = (i + 1)%cx;
		pj = (j + 1)%cy;
		pk = (k + 1)%cz;
		#Cell indexes to the left
		ni = (i - 1 + cx)%cx;
		nj = (j - 1 + cy)%cy;
		nk = (k - 1 + cz)%cz;
		#Neighbour cells
		cell = 24*(i + cx*(j + cy*k));
		right = 24*(pi + cx*(j + cy*k));
		left = 24*(ni + cx*(j + cy*k));
		front = 24*(i + cx*(pj + cy*k));
		back = 24*(i + cx*(nj + cy*k));
		up = 24*(i + cx*(j + cy*pk));
		down = 24*(i + cx*(j + cy*nk));
		diagp = 24*(pi + cx*(pj + cy*pk));
		diagn = 24*(ni + cx*(nj + cy*nk));
		#Unit cell list
		cellList = [];
		for p in range(0, 24):
			cellList.append([]);
		#end
		#Particle 0
		ID = 0;
		t = typeSi;
		r = [0, 0, 0];
		v = randomVelocity(sigmaSi);
		bonds = [down + 20, left + 5, back + 2, diagn + 23];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 1
		ID = 1;
		t = typeSi;
		r = [4, 4, 0];
		v = randomVelocity(sigmaSi);
		bonds = [cell + 3, cell + 4, down + 21, down + 22];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 2
		ID = 2;
		t = typeO;
		r = [1, 7, 1];
		v = randomVelocity(sigmaO);
		bonds = [cell + 6, front];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 3
		ID = 3;
		t = typeO;
		r = [3, 5, 1];
		v = randomVelocity(sigmaO);
		bonds = [cell + 1, cell + 6];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 4
		ID = 4;
		t = typeO;
		r = [5, 3, 1];
		v = randomVelocity(sigmaO);
		bonds = [cell + 1, cell + 7];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 5
		ID = 5;
		t = typeO;
		r = [7, 1, 1];
		v = randomVelocity(sigmaO);
		bonds = [cell + 7, right];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 6
		ID = 6;
		t = typeSi;
		r = [2, 6, 2];
		v = randomVelocity(sigmaSi);
		bonds = [cell + 2, cell + 3, cell + 8, cell + 9];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 7
		ID = 7;
		t = typeSi;
		r = [6, 2, 2];
		v = randomVelocity(sigmaSi);
		bonds = [cell + 4, cell + 5, cell + 10, cell + 11];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 8
		ID = 8;
		t = typeO;
		r = [1, 5, 3];
		v = randomVelocity(sigmaO);
		bonds = [cell + 6, cell + 12];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 9
		ID = 9;
		t = typeO;
		r = [3, 7, 3];
		v = randomVelocity(sigmaO);
		bonds = [cell + 6, front + 13];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 10
		ID = 10;
		t = typeO;
		r = [5, 1, 3];
		v = randomVelocity(sigmaO);
		bonds = [cell + 7, cell + 13];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 11
		ID = 11;
		t = typeO;
		r = [7, 3, 3];
		v = randomVelocity(sigmaO);
		bonds = [cell + 7, right + 12];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 12
		ID = 12;
		t = typeSi;
		r = [0, 4, 4];
		v = randomVelocity(sigmaSi);
		bonds = [cell + 8, cell + 14, left + 11, left + 17];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 13
		ID = 13;
		t = typeSi;
		r = [4, 0, 4];
		v = randomVelocity(sigmaSi);
		bonds = [cell + 10, cell + 15, back + 9, back + 16];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 14
		ID = 14;
		t = typeO;
		r = [1, 3, 5];
		v = randomVelocity(sigmaO);
		bonds = [cell + 12, cell + 18];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 15
		ID = 15;
		t = typeO;
		r = [3, 1, 5];
		v = randomVelocity(sigmaO);
		bonds = [cell + 13, cell + 18];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 16
		ID = 16;
		t = typeO;
		r = [5, 7, 5];
		v = randomVelocity(sigmaO);
		bonds = [cell + 19, front + 13];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 17
		ID = 17;
		t = typeO;
		r = [7, 5, 5];
		v = randomVelocity(sigmaO);
		bonds = [cell + 19, right + 12];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 18
		ID = 18;
		t = typeSi;
		r = [2, 2, 6];
		v = randomVelocity(sigmaSi);
		bonds = [cell + 14, cell + 15, cell + 20, cell + 21];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 19
		ID = 19;
		t = typeSi;
		r = [6, 6, 6];
		v = randomVelocity(sigmaSi);
		bonds = [cell + 16, cell + 17, cell + 22, cell + 23];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 20
		ID = 20;
		t = typeO;
		r = [1, 1, 7];
		v = randomVelocity(sigmaO);
		bonds = [cell + 18, up];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 21
		ID = 21;
		t = typeO;
		r = [3, 3, 7];
		v = randomVelocity(sigmaO);
		bonds = [cell + 18, up + 1];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 22
		ID = 22;
		t = typeO;
		r = [5, 5, 7];
		v = randomVelocity(sigmaO);
		bonds = [cell + 19, up + 1];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		#Particle 23
		ID = 23;
		t = typeO;
		r = [7, 7, 7];
		v = randomVelocity(sigmaO);
		bonds = [cell + 19, diagp];
		cellList[ID] = [t, r, v, cell + ID, bonds];
		return cellList;
	#end
	
	def writeUnitCell(i, j, k, ifile):
		R = b*np.array([i, j, k]);
		db = 0.125*b;
		cellList = unitCell(i, j, k);
		for l in cellList:
			#Type
			ifile.write(str(l[0]) + " ");
			#Position
			for d in range(0, 3):
				ifile.write(str(R[d] + db[d]*l[1][d]) + " ");
			#end
			#Velocity
			for d in range(0, 3):
				ifile.write(str(l[2][d]) + " ");
			#end
			#ID
			ifile.write(str(l[3]));
			# #bonds
			# for bond in l[4]:
			# 	ifile.write(" " + str(bond));
			# #end
			ifile.write("\n");
		#end
	#end

	#------#
	# main #
	#------#
	createDir(path);
	foldername = path + "/init";
	createDir(foldername);
	createDir(path + "/00000");
	initfilename = path + "/00000/combined.xyz";

	initfile = open(initfilename, "w");
	initfile.write(str(N) + "\ntype x y z vx vy vz ID {bonds}\n");

	for i in range(0, cx):
		for j in range(0, cy):
			for k in range(0, cz):
				writeUnitCell(i, j, k, initfile);
			#end
		#end
	#end
	initfile.close();
	removeVelocity(initfilename);
	return N;
#end

def argon(c, b, path, T, typeAr = 0):
	if(len(c) != 3):
		return;
	#end
	c = np.array(c);
	p = np.array(p);
	s = c*b;
	#Simulation size
	sx = s[0];
	sy = s[1];
	sz = s[2];
	#number of atoms in each direction
	cx = c[0];
	cy = c[1];
	cz = c[2];
	rcut = 3.0;
	numtypes = 1;
	#size of simulation
	s = b*np.array(c);
	#create directory
	createDir(path);
	foldername = path + "/init";
	createDir(foldername);
	initfilename = foldername + "/init.xyz";
	datafilename = foldername + "/data.dat";

	#Particle properties
	m = 39.948;
	q = 0;
	#LJ parameters
	R = 3.405;
	D = 0.010324;
	
	#standard deviation for velocity distribution
	sigma = 0.000911837*np.sqrt(T/m);

	def randomV():
		return random.gauss(0, sigma);
	#end
	
	def randomVelocity():
		return np.array([randomV(), randomV(), randomV()]);
	#end

	def saveToFile(r, v, ID, ifile):
		ifile.write(str(typeAr));
		ifile.write(" ");
		ifile.write(str(r[0]));
		ifile.write(" ");
		ifile.write(str(r[1]));
		ifile.write(" ");
		ifile.write(str(r[2]));
		ifile.write(" ");
		ifile.write(str(v[0]));
		ifile.write(" ");
		ifile.write(str(v[1]));
		ifile.write(" ");
		ifile.write(str(v[2]));
		ifile.write(" ");
		ifile.write(str(ID));
		ifile.write("\n");
	#end

	def unitcell(i, j, k, ifile):
		cell = i + cx*(j + cy*k);
		r0 = np.array([i*b, j*b, k*b]);
		r1 = r0 + np.array([0.5*b, 0.5*b, 0]);
		r2 = r0 + np.array([0.5*b, 0, 0.5*b]);
		r3 = r0 + np.array([0, 0.5*b, 0.5*b]);
		v0 = randomVelocity();
		v1 = randomVelocity();
		v2 = randomVelocity();
		v3 = randomVelocity();
		saveToFile(r0, v0, 4*cell, ifile);
		saveToFile(r1, v1, 4*cell + 1, ifile);
		saveToFile(r2, v2, 4*cell + 2, ifile);
		saveToFile(r3, v3, 4*cell + 3, ifile);
	#end

	#Make init file
	initfile = open(initfilename, "w");
	#Writes number of particles
	initfile.write(str(4*cx*cy*cz));
	initfile.write("\n<comment>\n");
	for i in range(0, cx):
		for j in range(0, cy):
			for k in range(0,  cz):
				unitcell(i, j, k, initfile);
			#end
		#end
	#end
	initfile.close();

	#Make data file
	#write datafile
	datafile = open(datafilename, "w");
	datafile.write("1 1 1\n");
	datafile.write(str(sx) + " " + str(sy) + " " + str(sz) + "\n");
	datafile.write(str(rcut*R) + "\n");
	datafile.write(str(numtypes) + "\n");
	#Particle properties
	datafile.write(str(typeAr) + " " + str(m) + " " + str(q) + "\n");
	datafile.write("#\n"); #LJ
	datafile.write(str(typeAr) + " " + str(R) + " " + str(D) + "\n");
	datafile.write("#\n"); #Radial
	datafile.write("#\n"); #Angular
	datafile.close();
#end



