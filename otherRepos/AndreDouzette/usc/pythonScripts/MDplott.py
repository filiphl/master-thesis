from __future__ import division;
import numpy as np;
import matplotlib.pyplot as plt;
import os;

def plotEnergy(path, V, K, t, label = ''):
	#Making plot directory
	dirname = runName + "/plot/";
	if(not os.path.exists(dirname)):
		os.makedirs(dirname);
	#end
	#Plotting kinetic energy
	plt.figure(100);
	if(label == ''):
		filename = path + "KineticEnergy.pdf";
	else:
		filename = path + "KineticEnergy_" + label + ".pdf";
	#end
	plt.plot(t, K, '-b');
	plt.title("Kinetic energy");
	plt.xlabel("Time[fs]");
	plt.ylabel("Energy[eV]");
	plt.savefig(filename, bbox_inches = 'tight');
	plt.close();
	#Plotting potential energy
	plt.figure(100);
	if(label == ''):
		filename = path + "PotentialEnergy.pdf";
	else:
		filename = path + "PotentialEnergy_" + label + ".pdf";
	#end
	plt.plot(t, V, '-g');
	plt.title("Potential energy");
	plt.xlabel("Time[fs]");
	plt.ylabel("Energy[eV]");
	plt.savefig(filename, bbox_inches = 'tight');
	plt.close();
	#Plotting total energy
	plt.figure(100);
	if(label == ''):
		filename = path + "TotalEnergy.pdf";
	else:
		filename = path + "TotalEnergy_" + label + ".pdf";
	#end
	plt.plot(t, V + K, '-r');
	error = np.std(V + K)/np.mean(V + K);
	plt.title("Total energy (Error: %.3g)"%error);
	plt.xlabel("Time[fs]");
	plt.ylabel("Energy[eV]");
	plt.savefig(filename, bbox_inches = 'tight');
	plt.close();
#end

def plotTemperature(path, T, t, label = ''):
	#Making plot directory
	dirname = runName + "/plot/";
	if(not os.path.exists(dirname)):
		os.makedirs(dirname);
	#end
	#Plotting temperature
	plt.figure(100);
	if(label == ''):
		filename = path + "Temperature.pdf";
	else:
		filename = path + "Temperature_" + label + ".pdf";
	#end
	plt.plot(t, T, '-b');
	plt.title("Temperature");
	plt.xlabel("Time[fs]");
	plt.ylabel("Temperature[K]");
	plt.savefig(filename, bbox_inches = 'tight');
	plt.close();
#end

def plotDisplacement(path, r2, t, label = '', indexes = [0, 1, 2]):
	#Making plot directory
	dirname = runName + "/plot/";
	if(not os.path.exists(dirname)):
		os.makedirs(dirname);
	#end
	indexes = np.array(indexes);
	#plotting <r^2>
	if(label == ''):
		filename = path + "Displacement.pdf";
	else:
		filename = path + "Displacement_" + label + ".pdf";
	#end
	plt.figure(100); plt.hold('on');
	for i in range(0, len(indexes)):
		plt.plot(t, r2[:, indexes[i]]);
	#end
	legend = ['Si', 'O', 'H'];
	plt.legend([legend[i] for i in indexes], loc = 2);
	plt.title(r"$<r^2>$");
	plt.xlabel("Time[fs]");
	plt.ylabel(r"Displacement$[\AA^2]$");
	plt.savefig(filename, bbox_inches = 'tight');
	plt.close();
#end

def plotRadialDistribution(path, g, r, label = '', indexes = [0, 1, 2]):
	#Making plot directory
	dirname = runName + "/plot/";
	if(not os.path.exists(dirname)):
		os.makedirs(dirname);
	#end
	indexes = np.array(indexes);
	plt.figure(100); plt.hold('on');
	for i in range(0, len(indexes)):
		plt.plot(t, r2[:, indexes[i]]);
	#end
	if(label == ''):
		filename = path + "RadialDistribution.pdf";
	else:
		filename = path + "RadialDistribution_" + label + ".pdf";
	#end
	legend = ['Si', 'O', 'H'];
	plt.legend([legend[i] for i in indexes], loc = 2);
	plt.title(r"$<r^2>$");
	plt.xlabel("Time[fs]");
	plt.ylabel(r"Displacement$[\AA^2]$");
	plt.savefig(filename, bbox_inches = 'tight');
	
	plt.figure(100);
	plt.hold('on');
	L = [];
	for i in range(0, len(indexes)):
		for j in range(0, i + 1):
			plt.plot(r, g[i, j]);
			L.append(legend[i] + "-" + legend[j]);
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
		for j in range(0, i + 1):
			if(label == ''):
				filename = path + "RadialDistribution" + legend[i] + legend[j] + ".pdf";
			else:
				filename = path + "RadialDistribution" + legend[i] + legend[j] + "_" + label + ".pdf";
			#end
			n += 1;
			plt.figure(100 + n);
			plt.plot(r, g[i, j]);
			plt.title("Radial distribution function (" + legend[i] + "-" + legend[j] + ")");
			plt.xlabel(r"Position[$\AA$]");
			plt.ylabel(r"Number density$[\AA^{-3}]$");
			plt.savefig(filename, bbox_inches = 'tight');
			plt.close();
		#end
	#end
#end