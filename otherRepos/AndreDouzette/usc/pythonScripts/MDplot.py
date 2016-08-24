from __future__ import division;
import numpy as np;
# import matplotlib.pyplot as plt;
noplt = False;
try:
	import matplotlib.pyplot as plt;
except ImportError:
	noplt = True;
#end
import os;


#Update stuff with noplt!

def saveData(filename, x, y):
	datafile = open(filename, 'w');
	#Write each row of the [x y z] matrix to each line in the file
	for i in range(0, len(x)):
		datafile.write("%g %g\n"%(x[i], y[i]));
	#end
	datafile.close();
#end

def plotEnergy(path, V, K, t, label = ''):
	#Making plot directory
	dirname = path + "/plot/";
	if(not os.path.exists(dirname)):
		os.makedirs(dirname);
	#end
	datadirname = path + "/data/";
	if(not os.path.exists(datadirname)):
		os.makedirs(datadirname);
	#end
	#Plotting kinetic energy
	if(label == ''):
		filename = path + "/plot/KineticEnergy.pdf";
		datafilename = path + "/data/KineticEnergy.dat";
	else:
		filename = path + "/plot/KineticEnergy_" + label + ".pdf";
		datafilename = path + "/data/KineticEnergy_" + label + ".dat";
	#end
	if(not noplt):
		plt.figure(100);
		plt.plot(t, K, '-b');
		plt.title("Kinetic energy");
		plt.xlabel("Time[fs]");
		plt.ylabel("Energy[eV]");
		plt.savefig(filename, bbox_inches = 'tight');
	#end
	saveData(datafilename, t, K);
	#Plotting potential energy
	if(label == ''):
		filename = path + "/plot/PotentialEnergy.pdf";
		datafilename = path + "/data/PotentialEnergy.dat";
	else:
		filename = path + "/plot/PotentialEnergy_" + label + ".pdf";
		datafilename = path + "/data/PotentialEnergy_" + label + ".dat";
	#end
	if(not noplt):
		plt.figure(101);
		plt.plot(t, V, '-g');
		plt.title("Potential energy");
		plt.xlabel("Time[fs]");
		plt.ylabel("Energy[eV]");
		plt.savefig(filename, bbox_inches = 'tight');
	#end
	saveData(datafilename, t, V);
	#Plotting total energy
	if(label == ''):
		filename = path + "/plot/TotalEnergy.pdf";
		datafilename = path + "/data/TotalEnergy.dat";
	else:
		filename = path + "/plot/TotalEnergy_" + label + ".pdf";
		datafilename = path + "/data/TotalEnergy_" + label + ".dat";
	#end
	error = np.std(V + K)/np.mean(V + K);
	if(not noplt):
		plt.figure(102);
		plt.plot(t, V + K, '-r');
		plt.title("Total energy (Error: %.3g)"%error);
		plt.xlabel("Time[fs]");
		plt.ylabel("Energy[eV]");
		plt.savefig(filename, bbox_inches = 'tight');
	#end
	saveData(datafilename, t, V + K);
	if(not noplt):
		plt.figure(100); plt.close();
		plt.figure(101); plt.close();
		plt.figure(102); plt.close();
	#end
#end

def plotTemperature(path, T, t, label = ''):
	#Making plot directory
	dirname = path + "/plot/";
	if(not os.path.exists(dirname)):
		os.makedirs(dirname);
	#end
	datadirname = path + "/data/";
	if(not os.path.exists(datadirname)):
		os.makedirs(datadirname);
	#end
	
	#Plotting temperature
	if(label == ''):
		filename = path + "/plot/Temperature.pdf";
		datafilename = path + "/data/Temperature.dat";
	else:
		filename = path + "/plot/Temperature_" + label + ".pdf";
		datafilename = path + "/data/Temperature_" + label + ".dat";
	#end
	if(not noplt):
		plt.figure(100);
		plt.plot(t, T, '-b');
		plt.title("Temperature");
		plt.xlabel("Time[fs]");
		plt.ylabel("Temperature[K]");
		plt.savefig(filename, bbox_inches = 'tight');
		plt.close();
	#end
	saveData(datafilename, t, T);
#end

def plotDisplacement(path, r2, t, label = ''):
	#Making plot directory
	dirname = path + "/plot/";
	if(not os.path.exists(dirname)):
		os.makedirs(dirname);
	#end
	datadirname = path + "/data/";
	if(not os.path.exists(datadirname)):
		os.makedirs(datadirname);
	#end
	legend = ['Si', 'O', 'H'];
	D = np.zeros(np.shape(r2));
	for i in range(0, 3):
		D[1:, i] = r2[1:, i]/(6*t[1:]);
	#end
	indexes = [];
	for i in range(0, 3):
		if(r2[:, i].any()):
			indexes.append(i);
		#end
	#end
	if(not indexes):
		#All arrays are zeros
		return;
	#end
	indexes = np.array(indexes);
	
	#plotting <r^2>
	if(label == ''):
		filename = path + "/plot/Displacement.pdf";
	else:
		filename = path + "/plot/Displacement_" + label + ".pdf";
	#end
	plt.figure(100); plt.hold('on');
	for i in range(0, len(indexes)):
		plt.plot(t, r2[:, indexes[i]]);
		if(label == ''):
			datafilename = path + "/data/Displacement_" + legend[indexes[i]] + ".dat";
		else:
			datafilename = path + "/data/Displacement_" + legend[indexes[i]] + "_" + label + ".dat";
		#end
		saveData(datafilename, t, r2[:, indexes[i]]);
	#end
	plt.legend([legend[i] for i in indexes], loc = 2);
	plt.title(r"$<r^2>$");
	plt.xlabel("Time[fs]");
	plt.ylabel(r"Displacement$[\AA^2]$");
	plt.savefig(filename, bbox_inches = 'tight');
	plt.close();
	
	#plotting diffusion
	if(label == ''):
		filename = path + "/plot/Diffusion.pdf";
	else:
		filename = path + "/plot/Diffusion_" + label + ".pdf";
	#end
	plt.figure(100); plt.hold('on');
	for i in range(0, len(indexes)):
		plt.plot(t, D[:, indexes[i]]);
		if(label == ''):
			datafilename = path + "/data/Diffusion_" + legend[indexes[i]] + ".dat";
		else:
			datafilename = path + "/data/Diffusion_" + legend[indexes[i]] + "_" + label + ".dat";
		#end
		saveData(datafilename, t, D[:, indexes[i]]);
	#end
	legend = ['Si', 'O', 'H'];
	plt.legend([legend[i] for i in indexes], loc = 1);
	plt.title("Diffusion coefficient");
	plt.xlabel("Time[fs]");
	plt.ylabel(r"$[\AA^2/fs]$");
	plt.savefig(filename, bbox_inches = 'tight');
	plt.close();
#end

def plotRadialDistribution(path, g, r, label = '', show = False):
	#Making plot directory
	dirname = path + "/plot/";
	if(not os.path.exists(dirname)):
		os.makedirs(dirname);
	#end
	datadirname = path + "/data/";
	if(not os.path.exists(datadirname)):
		os.makedirs(datadirname);
	#end
	
	#Finding which particle types are present
	indexes = [];
	for i in range(0, 3):
		for j in range(0, 3):
			if(g[i, j, :].any()):
				indexes.append(i);
				break;
			#end
		#end
	#end
	# print indexes;
	if(not indexes):
		#All arrays are zeros
		return;
	#end
	if(label == ''):
		filename = path + "/plot/RadialDistribution.pdf";
	else:
		filename = path + "/plot/RadialDistribution_" + label + ".pdf";
	#end
	
	legend = ['Si', 'O', 'H'];
	plt.figure(100);
	plt.hold('on');
	L = [];
	for i in range(0, len(indexes)):
		ii = indexes[i];
		for j in range(0, i + 1):
			jj = indexes[j];
			plt.plot(r, g[jj, ii]);
			L.append(legend[jj] + "-" + legend[ii]);
		#end
	#end
	plt.legend(L, loc = 1);
	plt.title("Radial distribution function");
	plt.xlabel(r"Position[$\AA$]");
	plt.ylabel(r"Number density$[\AA^{-3}]$");
	plt.savefig(filename, bbox_inches = 'tight');
	plt.close();
	
	n = 0;
	for i in range(0, len(indexes)):
		ii = indexes[i];
		for j in range(0, i + 1):
			jj = indexes[j];
			if(label == ''):
				filename = path + "/plot/RadialDistribution_" + legend[jj] + "-" + legend[ii] + ".pdf";
				datafilename = path + "/data/RadialDistribution_" + legend[jj] + "-" + legend[ii] + ".dat";
			else:
				filename = path + "/plot/RadialDistribution_" + legend[jj] + "-" + legend[ii] + "_" + label + ".pdf";
				datafilename = path + "/data/RadialDistribution_" + legend[jj] + "-" + legend[ii] + "_" + label + ".dat";
			#end
			n += 1;
			plt.figure(100 + n);
			plt.plot(r, g[jj, ii]);
			plt.title("Radial distribution function (" + legend[jj] + "-" + legend[ii] + ")");
			plt.xlabel(r"Position[$\AA$]");
			plt.ylabel(r"Number density$[\AA^{-3}]$");
			plt.savefig(filename, bbox_inches = 'tight');
			plt.close();
			saveData(datafilename, r, g[jj, ii]);
		#end
	#end
#end

def plotAngularDistribution(path, angleDistribution, angle, label = '', show = False):
	#Making plot directory
	dirname = path + "/plot/";
	if(not os.path.exists(dirname)):
		os.makedirs(dirname);
	#end
	datadirname = path + "/data/";
	if(not os.path.exists(datadirname)):
		os.makedirs(datadirname);
	#end
	
	legend = ['O-Si-O', 'Si-O-Si', 'Si-O-H', 'H-O-H'];
	indexes = [];
	for i in range(0, 4):
		if(angleDistribution[i, :].any()):
			indexes.append(i);
		#end
	#end
	if(not indexes):
		#All arrays are zeros
		return;
	#end
	#plotting anglular distribution
	if(label == ''):
		filename = path + "/plot/AngularDistribution.pdf";
	else:
		filename = path + "/plot/AngularDistribution_" + label + ".pdf";
	#end
	plt.figure(100); plt.hold('on');
	for ind in indexes:
		plt.plot(angle, angleDistribution[ind, :]);
	#end
	plt.legend([legend[i] for i in indexes], loc = 2);
	plt.title("Angle distribution");
	plt.xlabel("Angle");
	plt.ylabel("Angular density");
	plt.savefig(filename, bbox_inches = 'tight');
	plt.close();
	#Plotting as separate plots
	for ind in indexes:
		if(label == ''):
			filename = path + "/plot/AngularDistribution_" + legend[ind] + ".pdf";
			datafilename = path + "/data/AngularDistribution_" + legend[ind] + ".dat";
			
		else:
			filename = path + "/plot/AngularDistribution_" + legend[ind] + "_" + label + ".pdf";
			datafilename = path + "/data/AngularDistribution_" + legend[ind] + "_" + label + ".dat";
		#end
		plt.figure(100);
		plt.plot(angle, angleDistribution[ind, :]);
		plt.title("Angle distribution (" + legend[ind] + ")");
		plt.xlabel("Angle[deg]");
		plt.ylabel("Angular density");
		plt.savefig(filename, bbox_inches = 'tight');
		plt.close();
		saveData(datafilename, angle, angleDistribution[ind, :]);
	#end
#end

def plotForceDistribution(path, r, data, label = ''):
	#IMPLEMENT THIS
	pass;
#end

def plotSilanol(path, N, t, A, label = ''):
	#Making plot directory
	dirname = path + "/plot/";
	if(not os.path.exists(dirname)):
		os.makedirs(dirname);
	#end
	datadirname = path + "/data/";
	if(not os.path.exists(datadirname)):
		os.makedirs(datadirname);
	#end
	
	data = N/A;
	
	print np.shape(t), np.shape(data);
	
	#Plotting temperature
	if(label == ''):
		filename = path + "/plot/Silanol.pdf";
		datafilename = path + "/data/Silanol.dat";
	else:
		filename = path + "/plot/Silanol_" + label + ".pdf";
		datafilename = path + "/data/Silanol_" + label + ".dat";
	#end
	if(not noplt):
		plt.figure(100);
		plt.plot(t, data, '-b');
		plt.xlabel("Time[fs]");
		plt.ylabel("Number density[nm^{-2}]");
		plt.savefig(filename, bbox_inches = 'tight');
		plt.close();
	#end
	saveData(datafilename, t, data);
#end