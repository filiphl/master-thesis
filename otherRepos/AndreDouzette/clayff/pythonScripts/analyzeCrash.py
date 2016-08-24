from __future__ import division;
import matplotlib.pyplot as plt;
import numpy as np;
import simpleHole as hole;
from MDmanager import MDmanager as mdm;
from combineFiles import combiner;
import os;
import shutil;
import sys;
import multiprocessing as mpi;
import datetime;
import MDfilehandling as mdf;
import MDanalysis as mda;

import MDmanager as mdm; mdm = mdm.MDmanager();

if(len(sys.argv) > 1):
	runName = sys.argv[1];
else:
	print "Specify path to the simulation";
#end

dirs = [f for f in os.listdir(runName) if os.path.isdir(os.path.join(runName, f))];
dirs.remove('init');
dirs.sort();
dirs = [os.path.join(runName, d) for d in dirs];
#Find number of processors
datafile = open(runName + "/init/data.dat", 'r');
p = [int(i) for i in datafile.readline().split(' ')];
np = p[0]*p[1]*p[2];

#Find number of states with all processor files and energy file
stableStates = -1;
for i in range(0, len(dirs)):
	n = 0;
	for f in os.listdir(dirs[i]):
		if(f.endswith(".xyz") and f != "combined.xyz" and f != "sorted.xyz"):
			n += 1;
		#end
	#end
	if(n == np):
		stableStates += 1;
	else:
		print str(i) + " is unstable";
		break;
	#end
#end

mdm.setDir(runName + "/analyzed");
mdm.loadSimulation(runName, stableStates);
mdm.makeVmdData();
# mdm.visualize();
# exit();
# N = mdm.numberOfSilanol(states = range(0, stableStates, 500));
# plt.plot(N);
# plt.show();

# mdm.analyze();
# print "Making vmd input files";

mdm.radialDistribution(1024, label = "amorph", states = range(0,stableStates));
mdm.angularDistribution(1024, label = "amorph", states = range(0,stableStates));
