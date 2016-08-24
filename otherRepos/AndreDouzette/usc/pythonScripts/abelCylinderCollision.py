from __future__ import division;
import numpy as np;
import MDmanager as mdm; mdm = mdm.MDmanager();
import sys;
import os;
import MDdata; usc = MDdata.usc();
import shutil;

dt = 0.25;
dtsio2 = 1.25;

tauTherm = 50*dt;
tauNH    = 25*dt;

sim  = 100;
step = 12500;

simSample = 100;
stepSample = 7500;

simTherm = 100;
simStab  = 100;
stepTherm = 100;
stepStab  = 100;

# p = np.zeros(3);
# p[0] = sys.argv[1];
# p[1] = sys.argv[2];
# p[2] = sys.argv[3];

runName = sys.argv[1];
if(os.path.exists(runName)):
	print "Deleting old run";
	shutil.rmtree(runName);
#end

mdm.setProcessors([3, 2, 2]);
mdm.setDir(runName);
# mdm.setExecutable(sys.argv[5]);

#Number of sio2 cells in yz direction
cy = 20;
cz = 4;
#Cell thickness of wall
cw = 8;
#Cell thickness of collision object and water reservoir
cc = 6;

#Separation length between wall and lens systems in collision direction
sep1 = 10;
#Separation length in oposite direction
sep2 = 4;

c1 = np.array([cw, cy, cz]);
c2 = np.array([cc, cy, cz]);

#unit cell sizes
bsio2 = usc.sio2DensityToCellSize(2.196);
bh2o = usc.h2oDensityToCellSize(1.1); #Adding 20% due to separation between objects

#Simulation size
s1 = bsio2*c1;
s2 = bsio2*c2;

D = s2[0];
B = s2[1];

ch = np.array([int(s/bh2o + 0.5) for s in s2]);
st = s1 + [s2[0] + sep1 + sep2, 0, 0];

rcut = 5.5;
h = 0.8*s2[1];
d = 3*rcut;
w = 10;

indent = 4;

time = sim*step*dt;
position = 0.5*s2 + [0.5*s2[0] - 0.5*(d + w), 0, 0];
v = (sep1 + indent)/time;

#Making a lens of amarphous sio2
mdm.setTemperature(300);
mdm.setVolume(s2);
print "Creating SiO2 collision object"
mdm.makeCristobalite(c2);
mdm.simulate(1, 0.01, 100);
#Selecting the whole box of sio2
sel = mdm.selectBox([-1, -1, -1], [1000, 1000, 1000]);
print "Heating wall to 4000K";
mdm.setTemperature(4000);
mdm.thermalize(simTherm, dtsio2, tauNH, 10, stepTherm);
mdm.noseHoover(simStab, dtsio2, tauNH, stepStab);
print "Cooling down wall to 300K";
mdm.setTemperature(300);
mdm.thermalize(simTherm, dtsio2, tauTherm, 10, stepTherm);
sio2 = mdm.noseHoover(simStab, dtsio2, tauNH, stepStab);

co = mdm.makeBentCylinder(position, d, h, w = w);

mdm.makeVmdData();
mdm.visualize();
exit();

#Adding water surrounding the lens
print "Creating water";
mdm.setTemperature(300);
mdm.makeWater(ch);
selWater = mdm.selectBox([-1, -1, -1], [1000, 1000, 1000]);
print "Thermalizing water";
mdm.thermalize(simTherm, dt, tauTherm, 10, stepTherm);
mdm.noseHoover(simStab, dt, tauNH, stepStab);
water = mdm.makeCylinder(position, d + 3, h + 3, w = w, invert = True);
co = mdm.combine(co, water);
mdm.freeze(sel, v = [0, 0, 0]);
mdm.thermalize(simTherm, dt, tauTherm, 10, stepTherm);
mdm.noseHoover(simStab, dt, tauNH, stepStab);

#Creating the sio2 wall
mdm.setTemperature(300);
mdm.setVolume(s1);
print "Creating SiO2 wall";
sio2 = mdm.makeCristobalite(c1);
mdm.thermalize(simTherm, dtsio2, tauTherm, 10, stepTherm);
print "Heating wall to 4000K";
mdm.setTemperature(4000);
mdm.thermalize(simTherm, dtsio2, tauNH, 10, stepTherm);
mdm.noseHoover(simStab, dtsio2, tauNH, stepStab);
print "Cooling down wall to 300K";
mdm.setTemperature(300);
mdm.thermalize(simTherm, dtsio2, tauTherm, 10, stepTherm);
sio2 = mdm.noseHoover(simStab, dtsio2, tauNH, stepStab);

#Selecting the last 15% for freezing.
selWall = mdm.selectBox([s1[0] + 1, -1, -1], [0.85*s1[0], B + 1, B + 1]);
#Monitoring 1nm into the medium of the surface.
selMonitor = mdm.selectBox([-1, -1, -1], [10, B + 1, B + 1]);

#Combining the two systems
mdm.setVolume(st);
sio2 = mdm.moveParticles([s2[0] + sep1, 0, 0], sio2);
mdm.combine(co, sio2);
#Freezing objects
print "Equilibrating water with sio2 surfaces"
mdm.freeze(sel, v = [0, 0, 0]);
mdm.freeze(selWall, v = [0, 0, 0]);
# Equalibrating surface of wall with water
mdm.noseHoover(simStab, dt, tauNH, stepStab);

#Simulating the system, change to NH?
print "Moving object towards wall";
mdm.resetFreeze();
mdm.freeze(sel, v = [v, 0, 0]);
mdm.freeze(selWall, v = [0, 0, 0]);
mdm.noseHoover(sim, dt, tauNH, step);
#Sampling statistics
print "Sampling statistics";
mdm.resetFreeze();
mdm.freeze(sel, v = [0, 0, 0]);
mdm.freeze(selWall, v = [0, 0, 0]);
mdm.monitor(selMonitor);
mdm.noseHoover(simSample, dt, tauNH, stepSample);