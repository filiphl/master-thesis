#ifndef POTENTIAL_H
#define POTENTIAL_H

using namespace std;

#include "Particle.h"
#include <math.h>

//External libraries
#include <mpi.h>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>

class Potential{
	public:
		Potential();
		~Potential();
		void createForceFactor(double dt);
		void createVelocityFactor(double dt, double zeta);
		double force2(Particle* p);
		double force2Monitoring(Particle* p);
		double force3(Particle* p);
		double boundaryForce3(Particle* p);
		double twoBodyForce(Particle *p, int bond);
		double twoBodyOxygenForce(Particle *p, int bond);
		double twoBodyForceMonitoring(Particle *p, int bond);
		double twoBodyOxygenForceMonitoring(Particle *p, int bond);
		double threeBodyForce(Particle *p, int bond1, int bond2);
		double numberInterpolate(Particle *p, Particle *pb, double r2);
		void updateV(Particle *p);
		void updateVnh1(Particle *p, double zeta);
		void updateVnh2(Particle *p);
		void updateR(Particle *p, double dt);
		void print();
		double L0, E0, m0, v0, t0, T0, q0;
		double *m;
		double rc, rc2;
		int types, typeSi, typeOs, typeH, typeOh, typeO;
		double **Z, *Rb, *Db, *rInner2, *rOuter2, *Dbi;
		double **H, **D, **W, **eta, **r1si, **r4si, **r0, **r02, **xi;
		double ***B, ***c0;
	protected:
		void Tabulate2BPotential(int type1, int type2);
		void writePotentialFile(string name, double* r2, double* f, double* V, int n);
		double tau, taui, pi;
		int tableSize;
		double tableIncrement, tableIncrementInv;
		double oneThird;
		double ***fTable, ***VTable, *r2Table;
		double *fFactor, *vFactor;
		double *rv, *rv1, *rv2;
};

#endif