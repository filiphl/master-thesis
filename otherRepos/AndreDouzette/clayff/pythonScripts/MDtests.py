from __future__ import division;
import numpy as np;
import matplotlib.pyplot as plt;
import MDmanager as mdm;
import MDdata as mdd; usc = mdd.usc();

def runPotentialTest(path):
	datafilename = path + "/init/data.dat";
	datafile = open(datafilename, "r");
	line = datafile.readline();
	datafile.close();
	split = line.split(' ');
	p = 1;
	for s in split:
		p = p*int(float(s));
	#end
	
	#!Change this!
	#Path of executable file
	filename = os.path.expanduser("~/usc/bin/potentialTest");
	length = str(length);
	timesteps = str(timesteps);
	dt = str(dt);
	p = str(p);
	T = str(T);
	tau = str(tau);
	os.system("mpirun -np " + p + " " + filename + " " + path);
#end

def compile():
	print "Compiling";
	#Compiling c++ source
	os.system("cd ../src");
	os.system("mkdir -p path");
	#Making binary executables
	os.system("make clean -C ../src");
	os.system("make all -C ../src");
#end

def makeMolecule(t, r, v, filename):
	savefile = open(filename, 'w');
	savefile.write("%d\n"%len(t));
	savefile.write('PotentialTest\n');
	for i in range(0, len(t)):
		x = str(r[i, 0]);
		y = str(r[i, 1]);
		z = str(r[i, 2]);
		vx = str(v[i, 0]);
		vy = str(v[i, 1]);
		vz = str(v[i, 2]);
		savefile.write(str(t[i]) + " " + x + " " + y + " " + z + " " + vx + " " + vy + " " + vz + " " str(i) + "\n");
	#end
	savefile.close();
#end

def getSiO2TestMolecule(i):
	t = np.array([0, 1, 1]);
	v = np.array([0, 0, 0]);
	d = 0.125*sqrt(3)*b;
	rc = usc.rc;
	N = 5
	R = np.zeros(N, 3, 3);
	eps = 0.5*d;
	
	#cenetered in cell
	r = np.zeros(3, 3);
	r[0, 0] = eps;
	r[0, 1] = eps;
	r[0, 2] = eps;
	r[1, :] = r[0, :];
	r[1, 0] += d
	r[2, :] = r[0, :];
	r[2, 1] += d
	R[0, :, :] = r;
	#2 inner cells
	r = np.zeros(3, 3);
	r[0, 0] = rc - 0.5*d;
	r[0, 1] = eps;
	r[0, 2] = eps;
	r[1, :] = r[0, :];
	r[1, 0] += d
	r[2, :] = r[0, :];
	r[2, 1] += d
	R[1, :, :] = r;
	#3 inner cells
	r = np.zeros(3, 3);
	r[0, 0] = rc - 0.5*d;
	r[0, 1] = rc - 0.5*d;
	r[0, 2] = eps;
	r[1, :] = r[0, :];
	r[1, 0] += d
	r[2, :] = r[0, :];
	r[2, 1] += d
	R[2, :, :] = r;
	#1 outer cell
	r = np.zeros(3, 3);
	r[0, 0] = eps;
	r[0, 1] = eps - d;
	r[0, 2] = eps;
	r[1, :] = r[0, :];
	r[1, 0] += d
	r[2, :] = r[0, :];
	r[2, 1] += d
	R[3, :, :] = r;
	#2 outer cells
	r = np.zeros(3, 3);
	r[0, 0] = eps;
	r[0, 1] = eps - d;
	r[0, 2] = eps;
	r[1, :] = r[0, :];
	r[1, 0] -= d
	r[2, :] = r[0, :];
	r[2, 1] += d
	R[4, :, :] = r;
	
	#Enforce periodic conditions
	R = (R + 4*rc)%(4*rc);
	return t, R[i, :, :], v;
#end

def makeInitSio2TestState(path, i, s):
	foldername = path + "/init";
	path = os.path.expanduser("~/PotentialTest");
	if(os.path.exists(path)):
		shutil.rmtree(path);
	#end
	os.makedirs(path);
	foldername = foldername + "/init";
	if(os.path.exists(foldername)):
		shutil.rmtree(foldername);
	#end
	os.makedirs(foldername);
	#Molecule file
	t, r, v = getSiO2TestMolecule(i);
	makeMolecule(t, r, v, foldername + "/init.xyz");
	#Datafile
	datafilename = foldername + "/data.dat";
	datafile = open(datafilename, "w");
	datafile.write("1 1 1\n");
	datafile.write(str(s[0]) + " " + str(s[1]) + " " + str(s[2]) + "\n");
	datafile.close();
	mdf.split(path, [2, 2, 2], s);
#end

def getParticleData(foldername):
	
#end

def sio2potentialTest():
	b = usc.sio2DensityToCellSize(2.196);
	cp = 4;
	c = np.ones(3)*cp;
	s = c*(rc + 1e-6);
	compile();
	
	N = 5;
	for i in range(0, N):
		makeInitSio2TestState(path, i, s);
	#end
	
	#create sio2 molecule in outfile
	#run cpp test
	
#end