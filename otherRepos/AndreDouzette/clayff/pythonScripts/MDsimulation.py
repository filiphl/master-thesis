from __future__ import division;
from subprocess import call;
import os;
	
#Short and simple python wrapper for the c++ simulation program
def simulate(length, dt, path, timesteps = 1):
	datafilename = path + "/init/data.dat";
	#!Add some test to see if datafile exists
	datafile = open(datafilename, 'r');
	#Read number of processors
	line = datafile.readline();
	datafile.close();
	split = line.split(' ');
	p = 1;
	for s in split:
		p = p*int(float(s));
	#end
	
	#!Change this!
	#Path of executable file
	filename = os.path.expanduser("~/clayff/bin/simulate");
	length = str(length);
	timesteps = str(timesteps);
	dt = str(dt);
	p = str(p);
	# print "mpirun -np " + p + " " + filename + " " + path + " " + length + " " + timesteps + " " + dt;
	os.system("mpirun -np " + p + " " + filename + " " + path + " " + length + " " + timesteps + " " + dt);
#end

def thermalize(length, dt, tau, T, avglength, path, timesteps = 1):
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
	filename = os.path.expanduser("~/clayff/bin/thermalize");
	length = str(length);
	timesteps = str(timesteps);
	dt = str(dt);
	p = str(p);
	tau = str(tau);
	T = str(T);
	avglength = str(avglength);
	# print "mpirun -np " + p + " " + filename + " " + path + " " + length + " " + timesteps + " " + dt + " " + tau + " " + T;
	os.system("mpirun -np " + p + " " + filename + " " + path + " " + length + " " + timesteps + " " + dt + " " + tau + " " + avglength + " " + T);
#end

def noseHoover(length, dt, tau, T, path, timesteps = 1):
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
	filename = os.path.expanduser("~/clayff/bin/noseHoover");
	length = str(length);
	timesteps = str(timesteps);
	dt = str(dt);
	p = str(p);
	T = str(T);
	tau = str(tau);
	# print "mpirun -np " + p + " " + filename + " " + path + " " + length + " " + timesteps + " " + dt + " " + tau + " " + T;
	os.system("mpirun -np " + p + " " + filename + " " + path + " " + length + " " + timesteps + " " + dt + " " + tau + " " + T);
#end