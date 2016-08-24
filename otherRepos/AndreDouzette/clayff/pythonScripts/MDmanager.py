from __future__ import division;
import numpy as np;
import matplotlib.pyplot as plt;
import shutil;
import os;
import MDfilehandling as mdf;
import MDsimulation as mds;
import MDanalysis as mda;
import MDinit as mdi;
import MDmanipulation as mdm;
import MDplot as mdplot;

#python md API class
class MDmanager:
	simulations = [];
	stateNumbers = [];
	selections = 0;
	lastSim = -1;
	lastSimPath = "";
	lastState = 0;
	path = "";
	splitted = [];
	
	p = 2*np.ones(3, dtype = int);
	s = np.zeros(3);
	T = -1;
	
	def validSetup(this, sim = True):
		if(not sim):
			if(this.T < 0):
				print "Invalid temparature";
				return False;
			#end
			if(this.s[0] <= 0 or this.s[1] <= 0 or this.s[2] <= 0):
				print "Invalid volume";
				return False;
			#end
		#end
		if(this.p[0] <= 1 or this.p[1] <= 1 or this.p[2] <= 1):
			print "Invalid number of processes. Minimum process size: 2x2x2";
			return False;
		#end
		return True;
	#end
	
	def makeInitState(this, path, state, dest):
		mdf.combine(path, state, states = "last");
		splitDir = path + "/%05d"%state;
		#Removes all files in directory before recreating
		shutil.rmtree(dest);
		os.makedirs(dest);
		mdf.split(path, this.p, this.s, dest, splitDir + "/combined.xyz");
		#Get number of particles from datafile
		datafile = open(path + "/init/data.dat");
		lines = [l for l in datafile];
		N = int(lines[4]);
		return N;
	#end
	
	def reset(this):
		this.simulations = [];
		this.stateNumbers = [];
		this.lastSim = -1;
		this.lastSimPath = "";
		this.lastState = 0;
		this.splitted = [];
	#end
	
	def compile(this):
		print "Compiling";
		#Compiling c++ source
		os.system("cd ../src");
		os.system("mkdir -p path");
		#Making binary executables
		os.system("make clean -C ../src -s");
		os.system("make all -C ../src -s");
		# os.system("make all -C ../src"); #For bugtesting
	#end
	
	#number of processors
	def setProcessors(this, p):
		if(len(p) != 3):
			print "Number of processors needs to be specified in each dimension";
			return;
		#end
		this.p = np.array(p, dtype = int);
	#end
	
	#Simulation space
	def setVolume(this, s):
		if(len(s) != 3):
			print "Volume needs 3 components";
			return;
		#end
		this.s = np.array(s);
	#end
	
	def setTemperature(this, T):
		#
		this.T = T;
	#end
	
	def loadState(this, runPath, state):
		split = runPath.split("_");
		this.lastSim = int(split[0]);
		this.lastSimPath = runPath;
		this.lastState = state;
	#end
	
	def setDir(this, path):
		#
		this.path = path;
	#end
	
	def resetRun(this, folder):
		if(folder != ""):
			if(os.path.isdir(folder)):
				shutil.rmtree(folder);
			#end
		#end
	#end
	
	def makeWater(this, c, border = 0.8):
		#Border is in percentage of free space around water molecules
		path = this.path;
		p = this.p;
		s = this.s;
		T = this.T;
		resetRun = this.resetRun;
		if(path == ""):
			print "A save directory is needed to create a state";
			print "Use setDir(<dirName>) to set a save directory";
			return;
		#end
		if(s[0] == 0):
			print "Specify simulation volume";
			return;
		#end
		if(len(c) != 3 or c[0] == 0):
			print "Specify number of cells";
			return;
		#end
		#Removing all data from previous simulations
		lastSim = len(this.simulations) - 1;
		simPath = path + "/" + str(lastSim + 1) + "_H2OInit";
		resetRun(simPath);
		os.makedirs(simPath);
		N = mdi.water(s, c, simPath, T, border);
		
		datafilename = simPath + "/init/data.dat";
		datafile = open(datafilename, "w");
		datafile.write("1 1 1\n");
		datafile.write(str(s[0]) + " " + str(s[1]) + " " + str(s[2]) + "\n");
		datafile.write("0\n0\n");
		datafile.write(str(int(N)) + "\n");
		datafile.close();
		
		this.simulations.append(simPath);
		this.stateNumbers.append(0);
		this.splitted.append(False);
		return lastSim + 1;
	#end
	
	def makeCristobalite(this, c):
		resetRun = this.resetRun;
		path = this.path;
		p = this.p;
		s = this.s;
		T = this.T;
		if(path == ""):
			print "A save directory is needed to create a state";
			print "Use setDir(<dirName>) to set a save directory";
			exit();
		#end
		if(s[0] == 0):
			print "Specify simulation volume";
			exit();
		#end
		if(len(c) != 3 or c[0] == 0):
			print "Specify number of cells";
			exit();
		#end
		lastSim = len(this.simulations) - 1;
		simPath = path + "/" + str(lastSim + 1) + "_SiO2Init";
		resetRun(simPath);
		os.makedirs(simPath);
		N = mdi.cristobalite(s, c, T, simPath);
		
		datafilename = simPath + "/init/data.dat";
		datafile = open(datafilename, "w");
		datafile.write("1 1 1\n");
		datafile.write(str(s[0]) + " " + str(s[1]) + " " + str(s[2]) + "\n");
		datafile.write("0\n0\n");
		datafile.write(str(N) + "\n");
		datafile.close();
		
		this.simulations.append(simPath);
		this.stateNumbers.append(0);
		this.splitted.append(False);
		return lastSim + 1;
	#end
	
	def simulate(this, length, dt, timesteps, combine = False, resume = -1):
		lastSimPath = this.simulations[resume];
		lastState = this.stateNumbers[resume];
		lastSim = len(this.simulations) - 1;
		path = this.path
		if(path == ""):
			print "A save directory is needed to run simulation";
			print "Use setDir(<dirName>) to set a save directory";
			return;
		#end
		#Creating path of new simulation
		newSimPath = path + "/" + str(lastSim + 1) + "_simulation";
		#Copy old datafile to new location
		if(lastSim == -1):
			print "A initial state is needed for simulation";
			return;
		#end
		
		#Creating new directories and data
		mdf.makeDirs(length, newSimPath, reset = True);
		N = this.makeInitState(lastSimPath, lastState, newSimPath + "/00000");
		os.makedirs(newSimPath + "/init");
		datafilename = newSimPath + "/init/data.dat"
		datafile = open(datafilename, 'w');
		datafile.write(str(this.p[0]) + " " + str(this.p[1]) + " " + str(this.p[2]) + "\n");
		datafile.write(str(this.s[0]) + " " + str(this.s[1]) + " " + str(this.s[2]) + "\n");
		datafile.write(str(dt) + "\n");
		datafile.write(str(timesteps) + "\n");
		datafile.write(str(int(N)) + "\n");
		datafile.close();
		#Copy all other files in init directory
		initfiles = os.listdir(lastSimPath + "/init");
		for f in initfiles:
			if(f != "data.dat"):
				shutil.copyfile(lastSimPath + "/init/" + f, newSimPath + "/init/" + f);
			#end
		#end
		
		mds.simulate(length, dt, newSimPath, timesteps = timesteps);
		
		if(combine):
			print "Combining files";
			mdf.combine(newSimPath, length);
			print "Sorting files";
			mdf.sortStates(newSimPath, length);
		#end
		this.simulations.append(newSimPath);
		this.stateNumbers.append(length);
		this.splitted.append(True);
		return lastSim + 1;
	#end
	
	def thermalize(this, length, dt, tau, avglength, timesteps, combine = False, resume = -1):
		#Berendsen thermostat
		lastSimPath = this.simulations[resume];
		lastState = this.stateNumbers[resume];
		lastSim = len(this.simulations) - 1;
		T = this.T;
		path = this.path
		if(path == ""):
			print "A save directory is needed to run simulation";
			print "Use setDir(<dirName>) to set a save directory";
			return;
		#end
		#Creating path of new simulation
		newSimPath = path + "/" + str(lastSim + 1) + "_thermalization";
		#Copy old datafile to new location
		if(lastSim == -1):
			print "A initial state is needed for simulation";
			return;
		#end
		
		#Creating new directories and copy data
		mdf.makeDirs(length, newSimPath, reset = True);
		N = this.makeInitState(lastSimPath, lastState, newSimPath + "/00000");
		os.makedirs(newSimPath + "/init");
		datafilename = newSimPath + "/init/data.dat"
		datafile = open(datafilename, 'w');
		datafile.write(str(this.p[0]) + " " + str(this.p[1]) + " " + str(this.p[2]) + "\n");
		datafile.write(str(this.s[0]) + " " + str(this.s[1]) + " " + str(this.s[2]) + "\n");
		datafile.write(str(dt) + "\n");
		datafile.write(str(timesteps) + "\n");
		datafile.write(str(int(N)) + "\n");
		datafile.close();
		#Copy all other files in init directory
		initfiles = os.listdir(lastSimPath + "/init");
		for f in initfiles:
			if(f != "data.dat"):
				shutil.copyfile(lastSimPath + "/init/" + f, newSimPath + "/init/" + f);
			#end
		#end
		
		mds.thermalize(length, dt, tau, T, avglength, newSimPath, timesteps = timesteps);
		
		if(combine):
			print "Combining files";
			mdf.combine(length, newSimPath);
			print "Sorting files";
			mdf.sortStates(length, newSimPath);
		#end
		this.simulations.append(newSimPath);
		this.stateNumbers.append(length);
		this.splitted.append(True);
		return lastSim + 1;
	#end
	
	def noseHoover(this, length, dt, tau, timesteps, combine = False, resume = -1):
		lastSimPath = this.simulations[resume];
		lastState = this.stateNumbers[resume];
		lastSim = len(this.simulations) - 1;
		path = this.path
		T = this.T;
		if(path == ""):
			print "A save directory is needed to run simulation";
			print "Use setDir(<dirName>) to set a save directory";
			return;
		#end
		#Creating path of new simulation
		newSimPath = path + "/" + str(lastSim + 1) + "_noseHoover";
		#Copy old datafile to new location
		if(lastSim == -1):
			print "A initial state is needed for simulation";
			return;
		#end
		
		#Creating new directories and copy data
		mdf.makeDirs(length, newSimPath, reset = True);
		N = this.makeInitState(lastSimPath, lastState, newSimPath + "/00000");
		os.makedirs(newSimPath + "/init");
		datafilename = newSimPath + "/init/data.dat"
		datafile = open(datafilename, 'w');
		datafile.write(str(this.p[0]) + " " + str(this.p[1]) + " " + str(this.p[2]) + "\n");
		datafile.write(str(this.s[0]) + " " + str(this.s[1]) + " " + str(this.s[2]) + "\n");
		datafile.write(str(dt) + "\n");
		datafile.write(str(timesteps) + "\n");
		datafile.write(str(int(N)) + "\n");
		datafile.close();
		#Copy all other files in init directory
		initfiles = os.listdir(lastSimPath + "/init");
		for f in initfiles:
			if(f != "data.dat"):
				shutil.copyfile(lastSimPath + "/init/" + f, newSimPath + "/init/" + f);
			#end
		#end
		
		mds.noseHoover(length, dt, tau, T, newSimPath, timesteps = timesteps);
		
		if(combine):
			print "Combining files";
			mdf.combine(length, newSimPath);
			print "Sorting files";
			mdf.sortStates(length, newSimPath);
		#end
		this.simulations.append(newSimPath);
		this.stateNumbers.append(length);
		this.splitted.append(True);
		return lastSim + 1;
	#end
	
	def analyze(this, simulation = -1, label = ''):
		simPath = this.simulations[simulation];
		length = this.stateNumbers[simulation];
		#Read dt and step from file
		datafile = open(simPath + "/init/data.dat");
		datafile.readline(); datafile.readline();
		dt = float(datafile.readline());
		step = int(datafile.readline());
		datafile.close();
		t = step*dt*np.arange(0, length + 1);
		K, V, T, r2 = mda.getAnalysis(simPath, length);
		mdplot.plotEnergy(this.path, V, K, t, label = label);
		mdplot.plotTemperature(this.path, T, t, label = label);
		mdplot.plotDisplacement(this.path, r2, t, label = label);
		D = np.zeros(np.shape(r2));
		for i in range(0, 3):
			D[1:, i] = r2[1:, i]/(6*t[1:]);
		#end
		return K, V, T, r2, D, t;
	#end
	
	def calculateDisplacement(this, simulation = -1):
		path = this.simulations[simulation];
		length = this.stateNumbers[simulation];
		return mda.calculateDisplacement(path, length);
	#end
	
	def makeVmdData(this, simulation = -1):
		print "Making vmd input files";
		path = this.simulations[simulation];
		length = this.stateNumbers[simulation];
		mdf.makeVmdData(path, length);
		mdf.combine(path, length);
	#end
	
	def loadSimulation(this, path, length = None):
		#Finding length by number of directories minus one (for the data file directory)
		if(length == None):
			length = len(os.listdir(path)) - 2;
		#end
		#assuming state is split
		this.simulations.append(path);
		this.stateNumbers.append(length);
		this.splitted.append(True);
	#end
	
	def visualize(this):
		#Visualize using vmd
		mda.vmdVisualize();
	#end
	
	def getState(this):
		#
		return len(this.simulations) - 1;
	#end
	
	def radialDistribution(this, bins, simulation = -1, states = "all", label = ''):
		print "Calculating radial distribution";
		simPath = this.simulations[simulation];
		length = this.stateNumbers[simulation];
		mdf.combine(simPath, length, );
		r, g = mda.calculateRadialDistribution(simPath, bins, length, states);
		mdplot.plotRadialDistribution(this.path, g, r, label = label);
		return r, g;
	#end
	
	def moveParticles(this, dr, simulation = -1):
		length = this.stateNumbers[simulation];
		lastSim = this.simulations[simulation];
		splitted = this.splitted[simulation];
		path = this.path;
		simNumber = len(this.splitted) - 1;
		if(splitted):
			mdf.combine(lastSim, length, states = "last");
		#end
		loadfile = lastSim + "/%05d/combined.xyz"%length;
		newPath = path + "/" + str(simNumber + 1) + "_displaced";
		#Creating new directories and data
		if(not os.path.exists(newPath + "/init")):
			os.makedirs(newPath + "/init");
		#end
		datafilename = newPath + "/init/data.dat"
		datafile = open(datafilename, 'w');
		datafile.write(str(this.p[0]) + " " + str(this.p[1]) + " " + str(this.p[2]) + "\n");
		datafile.write(str(this.s[0]) + " " + str(this.s[1]) + " " + str(this.s[2]) + "\n");
		datafile.write("0\n0\n");
		#Number of particles
		lastDatafile = open(lastSim + "/init/data.dat");
		lines = [l for l in lastDatafile];
		N = int(lines[4]);
		datafile.write(str(N) + "\n");
		lastDatafile.close();
		datafile.close();
		#Copying mark files
		initfiles = os.listdir(lastSim + "/init");
		for f in initfiles:
			if(f != "data.dat"):
				shutil.copyfile(lastSimPath + "/init/" + f, newSimPath + "/init/" + f);
			#end
		#end
		#Creating directory to hold xyz files
		directory = newPath + "/00000";
		this.resetRun(directory);
		os.makedirs(directory);
		savefile = directory + "/combined.xyz";
		mdm.moveParticles(loadfile, savefile, dr, this.s);
		mdf.split(newPath, this.p, this.s, source = newPath + "/00000/combined.xyz");
		this.simulations.append(newPath);
		this.stateNumbers.append(0);
		this.splitted.append(True);
		return simNumber + 1;
	#end
	
	def displaceRadially(this, R, direction, simulation = -1):
		length = this.stateNumbers[simulation];
		lastSim = this.simulations[simulation];
		splitted = this.splitted[simulation];
		path = this.path;
		simNumber = len(this.splitted) - 1;
		if(splitted):
			mdf.combine(lastSim, length, states = "last");
		#end
		loadfile = lastSim + "/%05d/combined.xyz"%length;
		newPath = path + "/" + str(simNumber + 1) + "_radiallyDisplaced";
		#Creating new directories and data
		this.resetRun(newPath);
		os.makedirs(newPath + "/00000");
		os.makedirs(newPath + "/init");
		#Creatign data file
		datafilename = newPath + "/init/data.dat"
		datafile = open(datafilename, 'w');
		datafile.write(str(this.p[0]) + " " + str(this.p[1]) + " " + str(this.p[2]) + "\n");
		datafile.write(str(this.s[0]) + " " + str(this.s[1]) + " " + str(this.s[2]) + "\n");
		datafile.write("0\n0\n");
		lastDatafile = open(lastSim + "/init/data.dat");
		lines = [l for l in lastDatafile];
		N = int(lines[4]);
		datafile.write(str(N) + "\n");
		lastDatafile.close();
		datafile.close();
		initfiles = os.listdir(lastSimPath + "/init");
		for f in initfiles:
			if(f != "data.dat"):
				shutil.copyfile(lastSimPath + "/init/" + f, newSimPath + "/init/" + f);
			#end
		#end
		#Displacing particles
		mdm.displaceRadially(lastSim, newPath, R, direction)
		mdf.split(newPath, this.p, this.s, source = newPath + "/00000/combined.xyz");
		this.simulations.append(newPath);
		this.stateNumbers.append(0);
		this.splitted.append(True);
		return simNumber + 1;
	#end
	
	def combine(this, simulation1, simulation2):
		path = this.path;
		simNumber = len(this.splitted) - 1;
		
		length1 = this.stateNumbers[simulation1];
		lastSim1 = this.simulations[simulation1];
		splitted1 = this.splitted[simulation1];
		length2 = this.stateNumbers[simulation2];
		lastSim2 = this.simulations[simulation2];
		splitted2 = this.splitted[simulation2];
		if(splitted1):
			mdf.combine(lastSim1, length1, states = "last");
		#end
		if(splitted2):
			mdf.combine(lastSim2, length2, states = "last");
		#end
		newPath = path + "/" + str(simNumber + 1) + "_combined";
		#Creating new directories and data
		this.resetRun(newPath);
		os.makedirs(newPath + "/init");
		
		#Creating datafile
		datafilename = newPath + "/init/data.dat"
		datafile = open(datafilename, 'w');
		datafile.write(str(this.p[0]) + " " + str(this.p[1]) + " " + str(this.p[2]) + "\n");
		datafile.write(str(this.s[0]) + " " + str(this.s[1]) + " " + str(this.s[2]) + "\n");
		datafile.write("0\n0\n");
		#Number of particles
		lastDatafile = open(lastSim1 + "/init/data.dat");
		lines = [l for l in lastDatafile];
		N1 = int(lines[4]);
		lastDatafile.close();
		lastDatafile = open(lastSim2 + "/init/data.dat");
		lines = [l for l in lastDatafile];
		N2 = int(lines[4]);
		lastDatafile.close();
		datafile.write(str(N1 + N2) + "\n");
		datafile.close();
		
		#Combining files and splitting xyz file
		mdm.combine(lastSim1, lastSim2, newPath);
		mdf.split(newPath, this.p, this.s, source = newPath + "/00000/combined.xyz");
		this.simulations.append(newPath);
		this.stateNumbers.append(0);
		this.splitted.append(True);
		return simNumber + 1;
	#end
	
	#Returns angle distribution as following structure [O-Si-O, Si-O-Si, Si-O-H, H-O-H]
	def angularDistribution(this, bins, states = "all", simulation = -1, label = ''):
		print "Calculating angular distribution";
		length = this.stateNumbers[simulation];
		runName = this.simulations[simulation];
		mdf.combine(runName, length);
		angle, angleDistribution = mda.angularDistribution(bins, runName, length, states);
		mdplot.plotAngularDistribution(this.path, angleDistribution, angle*180/np.pi, label = label);
		return angle, angleDistribution;
	#end
	
	def makeLens(this, r0, d, h, direction = [1, 0, 0], invert = False, simulation = -1):
		length = this.stateNumbers[simulation];
		lastSim = this.simulations[simulation];
		splitted = this.splitted[simulation];
		path = this.path;
		simNumber = len(this.splitted) - 1;
		if(splitted):
			mdf.combine(lastSim, length, states = "last");
		#end
		loadfile = lastSim + "/%05d/combined.xyz"%length;
		newPath = path + "/" + str(simNumber + 1) + "_lens";
		#Creating new directories and data
		this.resetRun(newPath);
		os.makedirs(newPath + "/init");
		#Create datafile
		datafilename = newPath + "/init/data.dat"
		datafile = open(datafilename, 'w');
		datafile.write(str(this.p[0]) + " " + str(this.p[1]) + " " + str(this.p[2]) + "\n");
		datafile.write(str(this.s[0]) + " " + str(this.s[1]) + " " + str(this.s[2]) + "\n");
		datafile.write("0\n0\n");
		#Creating directory to hold xyz files
		N = mdm.makeLens(lastSim, newPath, r0, h, d, this.s, direction, invert);
		datafile.write(str(N) + "\n");
		datafile.close();
		mdf.split(newPath, this.p, this.s, source = newPath + "/00000/combined.xyz");
		this.simulations.append(newPath);
		this.stateNumbers.append(0);
		this.splitted.append(True);
		return simNumber + 1;
	#end
	
	def makeBentLens(this, r0, d, h, direction = [1, 0, 0], simulation = -1):
		R = (h*h + d*d)/(4*d);
		this.displaceRadially(R, -np.array(direction), simulation);
		return this.makeLens(r0, d, h, direction, invert = False, simulation = simulation);
	#end
	
	def toggleFreeze(this, simulation = -1, toggle = False):
		#Rename freeze file to unmark
		lastSim = this.simulations[simulation];
		f = lastSim + "/init/mark.dat";
		u = lastSim + "/init/unmark.dat";
		if(toggle):
			if(os.path.exists(u)):
				os.rename(u, f);
			#end
		else:
			if(os.path.exists(f)):
				os.rename(f, u);
			#end
		#end
	#end
	
	def numberOfSilanol(this, simulation = -1, states = "all"):
		length = this.stateNumbers[simulation];
		runName = this.simulations[simulation];
		return mda.silanolGroups(runName, length, states);
	#end
	
	def selectBox(this, r0, r1, simulation = -1, types = "all"):
		splitted = this.splitted[simulation];
		length = this.stateNumbers[simulation];
		lastSim = this.simulations[simulation];
		if(splitted):
			mdf.combine(lastSim, length, states = "last");
		#end
		markFile = lastSim + "/init/mark_%d.dat"%this.selections;
		this.selections += 1;
		filename = lastSim + "/%05d/combined.xyz"%length;
		mdm.selectBox(filename, markFile, r0, r1, types);
		return this.selections - 1;
	#end
	
	def monitor(this, selection = -1, simulation = -1):
		if(selection == -1):
			selection = this.selections - 1;
		#end
		lastSim = this.simulations[simulation];
		markFilename = lastSim + "/init/mark_%d.dat"%selection;
		if(not os.path.exists(markFilename)):
			print "%d is not a selection of the state %s"%(selection, lastSim);
			exit();
		#end
		monitoredFilename = lastSim + "/init/monitored.dat";
		shutil.copyfile(markFilename, monitoredFilename);
	#end
	
	def freeze(this, selection = -1, v = [0, 0, 0], simulation = -1):
		if(selection == -1):
			selection = this.selections - 1;
		#end
		if(selection == -1):
			print "There are no selections!";
			exit();
		#end
		lastSim = this.simulations[simulation];
		markFilename = lastSim + "/init/mark_%d.dat"%selection;
		if(not os.path.exists(markFilename)):
			print "%d is not a selection of the state %s"%(selection, lastSim);
			exit();
		#end
		frozenFilename = lastSim + "/init/frozen.dat";
		mdm.freeze(markFilename, frozenFilename, v);
	#end
	
	def resetFreeze(this, simulation = -1):
		lastSim = this.simulations[simulation];
		os.remove(lastSim + "/init/frozen.dat");
	#end
	
	def forceDistribution(this, bins, senderMark, direction, simulation = -1, states = "all", label = ''):
		length = this.stateNumbers[simulation];
		runName = this.simulations[simulation];
		r, data = mda.forceDistribution(runName, bins, senderMark, states, direction, length);
		fz = data[0, :];
		fr = data[1, :];
		ft = data[2, :];
		plt.figure(10);
		plt.plot(r, fz, '-b');
		plt.figure(11);
		plt.plot(r, fr, '-g');
		plt.figure(12);
		plt.plot(r, ft, '-r');
		plt.show();
	#end
#end










