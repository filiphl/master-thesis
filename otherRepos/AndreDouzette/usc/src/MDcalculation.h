#ifndef MDCALCULATION_H
#define MDCALCULATION_H

//include classes
#include "Particle.h"
#include "ParticleList.h"
#include "ParticleContainer.h"
#include "Potential.h"
//External libraries
#include <mpi.h>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
#include <string.h>

// using std::vector;
using namespace std;

class MDcalculation{
	public:
		MDcalculation(string path, int loadstate);
		void run(int Nt, double dt, string runName, int saveInterval);
		void thermalize(int Nt, double dt, double tau, int avglength, double T, string runName, int saveInterval);
		void noseHoover(int Nt, double dt, double Q, double T, string runName, int saveInterval, int Nchain);
		void readDatafile(string filename);
		void finalize();
		void readFile(string filename, int n);
		void saveFile(string filename, int n);
	private:
		void saveInteratomicForces(string runName, int n);
		void readMonitoredFile(string name);
		void addOxygenNeighbour(Particle* p, Particle* q);
		void addNeighbour(Particle* p, Particle* q);
		void addCloseBond(Particle* p, Particle* q);
		void findNeighbours();
		void readMarkedFile(string filename);
		void deleteLists();
		int coordsToRank(int x, int y, int z);
		void interchangeParticles();
		void interchangeParticles(int dim);
		void shiftCoordinatesOut(double* data, int nData, int dim);
		void shiftCoordinatesIn(double* data, int nData, int dim);
		void shiftBoundaryCoordinatesOut(double* data, int nData, int dim);
		void shiftBoundaryCoordinatesIn(double* data, int nData, int dim);
		void communicateBorder();
		void exportBoundary(ParticleList* rightExport, ParticleList* leftExport, int dim);
		void addBoundary(int* stats, double* data, int nStasts, int nData);
		void addToBorder(Particle* p);
		void addParticles(int* stats, double* data, int nStats, int nData);
		void addToSystem(Particle* p);
		void swapParticles();
		void integrate(double dt, bool monitoring);
		void integrateNoseHoover(double dt, double T, bool monitoring);
		void calculateForces(bool monitoring);
		void externalForce(Particle* p);
		void calculate2BodyForces(int i, int j, int k);
		void stringToArray(string str, string** strarray, int* length);
		void cellReact(ParticleList* list1, ParticleList* list2);
		void selfReact(ParticleList* list);
		void radialForce(Particle* p1, Particle* p2);
		void resetAllBonds();
		void angleForce(Particle* p1, Particle* pm, Particle* p2);
		void thermostat(double dt, double T0, double tau, double T);
		void calculateKineticEnergy();
		void calculateMeanDisplacement();
		bool bonded(Particle* p1, Particle* p2);
		void saveEnergy(string name, int n);
		Particle* searchForParticle(int id, ParticleList* list);
		Particle* findInsideParticle(int id, int i, int j, int k);
		Particle* findBoundaryParticle(int id, int i, int j, int k);
		Particle* findParticle(int id, int i, int j, int k);
		void boundary3BodyForce(int i, int j, int k);
		void calculate3BodyForces(int i, int j, int k);
		void calculate1BodyForces(int i, int j, int k);
		bool shareBondAngle(Particle* p1, Particle* p2);
		void boundary2Body();
		void createDataStructures();
		void convertTime(double seconds, char* timeleft);
		void createNewBonds();
		void addOxygenToList(ParticleList* list, int i, int j, int k);
		//Testing
		int globalParticles();
		int localParticles();
		void communicateWeights();
		void exportWeights(ParticleList* rightExport, ParticleList* leftExport, int dim);
		void addWeights(int* stats, double* data, int nStats, int nData);
		void calculateWeights();
		double interpolate(Particle* p, Particle* pb);
		void oxygenInteraction(Particle* p1, Particle* p2);
		void radialForceOxygen(Particle* p1, Particle* p2, double f);
		void sortParticles(Particle** data, int length);
		void quickSort(Particle** data, int left, int right);
		void findNeighboursInternal(int i0, int j0, int k0, int i1, int j1, int k1);
		void findNeighboursBetween(int i0, int j0, int k0, int i1, int j1, int k1);
		void findNeighboursBoundary(int i0, int j0, int k0, int i1, int j1, int k1);
		void virial(double y1, double y2, double w);
		void calculatePressure();
	private:
		//Processor rank
		int rank;
		//Number of processors
		int procs;
		//processor ids
		int xp, xn;
		int yp, yn;
		int zp, zn;
		//Number of cell lists in each dimension
		int Nx, Ny, Nz;
		//Size of unit cells
		double Lx, Ly, Lz;
		//Global simulation space
		double sx, sy, sz;
		//x, y, z processor coordinates
		int xrank, yrank, zrank;
		//number of processors in each dimension
		int procsx, procsy, procsz;
		//cell lists
		ParticleList**** cells;
		//Export lists
		ParticleList** exportList;
		//Communication lists
		ParticleList* communicateRight;
		ParticleList* communicateLeft;
		//Neighbouring processors
		int* neighbour;
		//Numerical conversion values
		double v0, L0, t0, T0;
		//Particle types
		int typeSi, typeO, typeH, typeOs, typeOh;
		//Force parameters
		//Particle properties
		//Radial
		double** H;
		double** D;
		double** W;
		double** nu;
		double** r1s;
		double** r4s;
		//Angular
		double*** tA;
		double*** kA;
		//Total number of particles
		int Nparticles;
		int* NparticlesType;
		//Energy
		double V;
		double K;
		//Displacement <r^2>
		double* rMeanSquared;
		double* rMeanSquaredSum;
		//Array of particle pointers
		ParticleContainer* container;
		//Potential and parameter container
		Potential* potential;
		ParticleList* moveList;
		double dr[3];
		double rc2;
		double rc;
		double timeForce, timeInt, timeCom, timeOther;
		double b3bf2;
		double* b3bf;
		bool monitoring;
		int NnotFrozen;
		
		vector<double> dot0;
		vector<double> dot1;
		vector<int> ind0;
		vector<int> ind1;
		double* zeta;
		int Nchain;
		double Q0, Q0i;
		double Qj, Qji;
		double *mA, *mB, *mC, *mD;
		double* pressure;
		int prN;
		double prL, prLi;
		int cy0;
		int pressureSamples;
		int globalprN;
		double pressure0;
};
#endif


















