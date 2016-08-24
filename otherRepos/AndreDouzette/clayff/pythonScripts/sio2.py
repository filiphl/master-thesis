from __future__ import division;
import numpy as np;
import simpleHole as hole;
from MDmanager import MDmanager as mdm;
from combineFiles import combiner;
import os;
import shutil;
import sys;
import multiprocessing as mpi;
import datetime;
import MDdata as mdd; usc = mdd.usc(); 

#Length (in timesteps) of simulation
sim = 1200;
step = 250;
#Timestep length (fs)
dt = 0.25;

therm = 100;
thermstep = 250;
stab = 300;
stabstep = 250

dateAndTime = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S");
runName = dateAndTime
if(len(sys.argv) > 1):
	runName = sys.argv[1];
#end
path = os.path.expanduser("~/clayff_sio2h2oSim/" + runName);
#Removes path if exists, resetting data
if(os.path.exists(path)):
	shutil.rmtree(path);
#end

#processors
if(mpi.cpu_count() == 8):
	#Cpu geometry
	p = [2, 2, 2];
elif(mpi.cpu_count() == 12):
	p = [3, 2, 2];
else:
	print "Unsupported number(" + str(mpi.cpu_count()) + ") of cpu cores, needs 8 or 12";
	exit();
#end

c = 8;

b = usc.sio2DensityToCellSize(2.196);
b = b;
c = np.array([c, c, c]);
s = b*c;

mdm = mdm();
mdm.compile();
mdm.setDir(path);
mdm.setProcessors(p);
mdm.setVolume(s);
mdm.setTemperature(4000);

mdm.makeCristobalite(c);

T = [4000, 3500, 3000, 2500, 2250, 2000, 1750, 1500, 1000, 500, 300];

for temp in T:
	print "T = ", temp, "K"
	mdm.setTemperature(temp);
	mdm.thermalize(therm, dt, 25*dt, 10, thermstep);
	mdm.makeVmdData();
	mdm.analyze(label = "%d_therm"%temp);
	mdm.noseHoover(stab, dt, 25*dt, stabstep);
	mdm.makeVmdData();
	mdm.analyze(label = "%d_stab"%temp);
#end
mdm.noseHoover(sim, dt, 25*dt, step);

# mdm.radialDistribution(1024, label = "amorph", states = range(sim-500,sim + 1));
# mdm.angularDistribution(1024, label = "amorph", states = range(sim-500,sim + 1));
# mdm.visualize();