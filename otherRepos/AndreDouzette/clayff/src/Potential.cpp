#include "Potential.h"
#include <iostream>
#include <mpi.h>
using namespace std;

//The CLAYFF model
Potential::Potential(){
	tau = 6.2831853071795864769252867665590057683943387987502116; //truncate
	taui = 1/tau;
	pi = 0.5*tau;
	
	//Numeric conversion factors. A = A0 An
	//These needs confirmation
	//Length: 1Aa
	L0 = 1;
	//Energy: 1eV
	E0 = 103.6382;
	//Mass: 1amu
	m0 = 1;
	//Velocity: sqrt(E0/m0) A/fs
	v0 = 1;
	//Time: L0/v0 fs
	t0 = 1;
	//Temperature: E0/kb K
	T0 = 1/(8.3148e-7);
	
	//Indexes types
	typeSi = 0;
	typeO = 1;
	typeH = 2;
	typeOs = 1;
	typeOh = 3;
	
	//Number of particle types
	types = 4;
	
	//2BF cutoff
	rc  = 5.5; //10
	
	double* z = new double[types];
	
	fFactor = new double[types];
	vFactor = new double[types];
	
	Z = new double*[types];
	m = new double[types];
	eps = new double*[types];
	sig = new double*[types];
	kt = new double**[types];
	theta0 = new double**[types];
	for(int i = 0; i < types; i++){
		Z[i] = new double[types];
		eps[i] = new double[types];
		sig[i] = new double[types];
		kt[i] = new double*[types];
		theta0[i] = new double*[types];
		m[i] = 0;
		z[i] = 0;
		for(int j = 0; j < types; j++){
			kt[i][j] = new double[types];
			theta0[i][j] = new double[types];
			eps[i][j] = 0;
			sig[i][j] = 0;
			for(int k = 0; k < types; k++){
				kt[i][j][k] = 0;
				theta0[i][j][k] = 0;
			}
		}
	}
	
	//Coulomb shielding
	double r1s = 5.649;
	double Dy = 1.07;
	r1si = 1/r1s;
	
	//Charge
	z[typeSi] = 2.1;
	z[typeOs] = -1.05;
	z[typeH]  = 0.32983;
	z[typeOh] = -0.65966;
	//Mass
	m[typeSi] = 28.0855;
	m[typeO]  = 15.9994;
	m[typeH]  = 1.00794;
	m[typeOh] = 15.9994;
	
	//Bond
	kr = 554.1349/(E0*L0);
	r0 = 1/L0;
	
	// double kcalpermol = 23.06;
	double kcalpermol = 0.0433634;
	
	//Angle
	kt[typeSi][typeOs][typeSi] = 0;
	kt[typeSi][typeOs][typeH]  = 30;
	kt[typeH][typeOh][typeH]   = 45.7696;
	
	theta0[typeSi][typeOs][typeSi] = 0;
	theta0[typeSi][typeOs][typeH]  = 109.47/180*pi;
	theta0[typeH][typeOh][typeH]   = 109.47/180*pi;
	
	kt[typeH][typeOs][typeSi] = kt[typeSi][typeOs][typeH];
	theta0[typeH][typeOs][typeSi] = theta0[typeSi][typeOs][typeH];
	
	//LJ
	eps[typeSi][typeSi] = 1.8405*0.000001;
	eps[typeOs][typeOs] = 0.1554;
	eps[typeH][typeH]   = 0;
	eps[typeOh][typeOh] = 0.1554;
	
	sig[typeSi][typeSi] = 3.7064;
	sig[typeOs][typeOs] = 3.5532;
	sig[typeH][typeH]   = 0;
	sig[typeOh][typeOh] = 3.5532;
	
	//Check this value!
	double z0 = 14.39964;
	
	for(int i = 0; i < types; i++){
		for(int j = 0; j < types; j++){
			eps[i][j] = eps[i][j]*kcalpermol;
			for(int k = 0; k < types; k++){
				kt[i][j][k] = kt[i][j][k]*kcalpermol;
			}
		}
	}
	for(int i = 0; i < types; i++){
		for(int j = 0; j < types; j++){
			Z[i][j] = z[i]*z[j]*z0;
			eps[i][j] = sqrt(eps[i][i]*eps[j][j]);
			sig[i][j] = 0.5*(sig[i][i] + sig[j][j]);
		}
	}
	
	for(int i = 0; i < types; i++){
		for(int j = 0; j < types; j++){
			Z[i][j] = Dy*Z[i][j]/E0;
			eps[i][j] = eps[i][j]/E0;
			sig[i][j] = sig[i][j]/L0;
			for(int k = 0; k < types; k++){
				kt[i][j][k] = kt[i][j][k]/E0;
			}
		}
	}
	
	rc = rc/L0;
	rc2 = rc*rc;
	delete[] z;
	
	//Calculating table elements
	tableSize = 1<<13;
	tableIncrement = rc2/(tableSize - 1);
	oneThird = 1.0/3.0;
	tableIncrementInv = 1/tableIncrement;
	fTable = new double**[types];
	VTable = new double**[types];
	r2Table = new double[tableSize];
	
	for(int i = 0; i < tableSize; i++){
		r2Table[i] = i*tableIncrement;
	}
	
	
	for(int i = 0; i < types; i++){
		fTable[i] = new double*[types];
		VTable[i] = new double*[types];
		for(int j = 0; j < types; j++){
			Tabulate2BPotential(i, j);
		}
	}
	
	rv = new double[3];
	rv1 = new double[3];
	rv2 = new double[3];
}

Potential::~Potential(){
	delete[] r2Table;
	delete[] m;
	delete[] fFactor;
	delete[] vFactor;
	//Add deletion for all arrays
	for(int i = 0; i < types; i++){
		for(int j = 0; j <= i; j++){
			delete[] fTable[i][j];
			delete[] VTable[i][j];
		}
		delete[] fTable[i];
		delete[] VTable[i];
	}
	delete[] fTable;
	delete[] VTable;
}

void Potential::Tabulate2BPotential(int type1, int type2){
	double r, r2, rs2, rs6, rs12;
	double f1, f2, f3, f4;
	double V1, V2, V3, V4;
	double Vc, fc;
	int n = tableSize;
	//Importing correct particle parameters
	// double r1si = this->r1si[type1][type2];
	double Z = this->Z[type1][type2];
	double eps = this->eps[type1][type2];
	double sig = this->sig[type1][type2];
	//Calculating potential in nodes given in rTable
	fTable[type1][type2] = new double[n];
	VTable[type1][type2] = new double[n];
	double* f = fTable[type1][type2];
	double* V = VTable[type1][type2];
	
	for(int i = 0; i < n; i++){
		r2 = r2Table[i];
		r = sqrt(r2);
		
		rs2 = sig*sig/r2;
		rs6 = rs2*rs2*rs2;
		rs12 = rs6*rs6;
		
		V[i] = Z*exp(-r*r1si)/r + eps*(rs12 - 2*rs6);
		f[i] = Z*exp(-r*r1si)/r2*(1/r + r1si) + 12*eps*(rs12 - rs6)/r2;
	}
	// return;
	//Truncating
	Vc = V[n - 1];
	fc = f[n - 1];
	for(int i = 0; i < n; i++){
		r = sqrt(r2Table[i]);
		V[i] = V[i] - Vc + (r - rc)*rc*fc; //f is stored as F/r
		f[i] = f[i] - fc*rc/r;
	}
}

double Potential::bondedForce(Particle *p, int bond){
	double r2, r, dr, f, Fv;
	Particle* q = p->bonds[bond];
	r2 = 0;
	for(int d = 0; d < 3; d++){
		rv[d] = q->r[d] - p->r[d];
		r2 += rv[d]*rv[d];
	}
	r = sqrt(r2);
	dr = r - r0;
	f = kr*dr/r;
	for(int d = 0; d < 3; d++){
		Fv = f*rv[d];
		p->F[d] +=  Fv;
		q->F[d] += -Fv;
	}
	return 0.5*kr*dr*dr;
}

//Calculates forces based on table lookup, and returns potential
double Potential::twoBodyForce(Particle *p, int bond){
	Particle* q = p->neighbours[bond];
	double t;
	int n;
	double *ft = fTable[p->type][q->type];
	double *Vt = VTable[p->type][q->type];
	double f, V, Fv;
	//Find interval in which to interpolate
	double r2 = 0;
	for(int d = 0; d < 3; d++){
		rv[d] = q->r[d] - p->r[d];
		r2 += rv[d]*rv[d];
	}
	
	r2 = tableIncrementInv*r2;
	n = (int)r2;
	//Interpolation scheme
	t = r2 - n;
	f = ft[n] + t*(ft[n + 1] - ft[n]);
	V = Vt[n] + t*(Vt[n + 1] - Vt[n]);
	//Covnerting to vectors and adding forces to particles(f = F/r)
	for(int d = 0; d < 3; d++){
		Fv = f*rv[d];
		p->F[d] += -Fv;
		q->F[d] +=  Fv;
	}
	
	if(q->outside){
		return 0.5*V;
	}
	return V;
}

double Potential::twoBodyForceMonitoring(Particle *p, int bond){
	Particle* q = p->neighbours[bond];
	double t;
	int n;
	double *ft = fTable[p->type][q->type];
	double *Vt = VTable[p->type][q->type];
	double f, V, Fv;
	//Find interval in which to interpolate
	double r2 = 0;
	for(int d = 0; d < 3; d++){
		rv[d] = q->r[d] - p->r[d];
		r2 += rv[d]*rv[d];
	}
	
	r2 = tableIncrementInv*r2;
	n = (int)r2;
	//Interpolation scheme
	t = r2 - n;
	f = ft[n] + t*(ft[n + 1] - ft[n]);
	V = Vt[n] + t*(Vt[n + 1] - Vt[n]);
	//Covnerting to vectors and adding forces to particles(f = |F|/r)
	if(p->monitored){
		p->index2B.push_back(q->id);
		if(q->monitored){
			q->index2B.push_back(p->id);
			for(int d = 0; d < 3; d++){
				Fv = f*rv[d];
				p->F[d] += -Fv;
				q->F[d] +=  Fv;
				p->forces2B.push_back(-Fv);
				q->forces2B.push_back(Fv);
			}
		}else{
			for(int d = 0; d < 3; d++){
				Fv = f*rv[d];
				p->F[d] += -Fv;
				q->F[d] +=  Fv;
				p->forces2B.push_back(-Fv);
			}
		}
	}else{
		if(q->monitored){
			q->index2B.push_back(p->id);
			for(int d = 0; d < 3; d++){
				Fv = f*rv[d];
				p->F[d] += -Fv;
				q->F[d] +=  Fv;
				q->forces2B.push_back(Fv);
			}
		}else{
			for(int d = 0; d < 3; d++){
				Fv = f*rv[d];
				p->F[d] += -Fv;
				q->F[d] +=  Fv;
			}
		}
	}
	
	// if(p->outside || q->outside){
	if(q->outside){
		return 0.5*V;
	}
	return V;
}

double Potential::threeBodyForce(Particle *p, int bond1, int bond2){
	Particle* q1 = p->bonds[bond1];
	Particle* q2 = p->bonds[bond2];
	
	//Retrieving parameters
	double t0 = this->theta0[q1->type][p->type][q2->type];
	double kt = this->kt[q1->type][p->type][q2->type];
	double s, f0;
	double Vt, a1, a2, b1, b2;
	double F1, F2;
	
	//Calculating cos of angle
	double r1 = 0;
	double r2 = 0;
	double c = 0;
	for(int d = 0; d < 3; d++){
		rv1[d] = q1->r[d] - p->r[d];
		rv2[d] = q2->r[d] - p->r[d];
		r1 += rv1[d]*rv1[d];
		r2 += rv2[d]*rv2[d];
		c += rv1[d]*rv2[d];
	}
	
	r1 = sqrt(r1);
	r2 = sqrt(r2);
	c = c/(r1*r2);
	
	if(c > 1){
		c = 1;
	}else if(c < -1){
		c = -1;
	}
	s = sqrt(1 - c*c);
	if(s < 0.001) s = 0.001;
	f0 = 2*kt*(acos(c) - t0)/s;
	
	//Force matrix components
	a1 = -f0*c/r1;
	b1 = f0/r2;
	a2 = f0/r1;
	b2 = -f0*c/r2;
	
	//Add collection of monitored particle forces?
	for(int d = 0; d < 3; d++){
		F1 = a1*rv1[d] + b1*rv2[d];
		F2 = a2*rv1[d] + b2*rv2[d];
		q1->F[d] += F1;
		q2->F[d] += F2;
		p->F[d] += -(F1 + F2);
	}
	//Potential energy calculation
	return f0*s*(1 - oneThird*(p->outside + q1->outside + q2->outside));
}

double Potential::force2(Particle* p){
	double V = 0;
	int N = p->neighbours.size();
	for(int i = 0; i < N; i++){
		V += twoBodyForce(p, i);
	}
	return V;
}

double Potential::force2Monitoring(Particle* p){
	double V = 0;
	int N = p->neighbours.size();
	for(int i = 0; i < N; i++){
		V += twoBodyForceMonitoring(p, i);
	}
	return V;
}

double Potential::force3(Particle* p){
	//Skipping particles without valance angles
	if(p->bonds.size() < 2){
		return 0;
	}
	double V = 0;
	int N = p->bonds.size();
	for(int i1 = 1; i1 < N; i1++){
		for(int i2 = 0; i2 < i1; i2++){
			V += bondedForce(p, i1);
			V += bondedForce(p, i2);
			V += threeBodyForce(p, i1, i2);
		}
	}
	return V;
}

double Potential::boundaryForce3(Particle* p){
	//Skipping particles without valance angles
	if(p->bonds.size() < 2){
		return 0;
	}
	double V = 0;
	int N = p->bonds.size();
	// cout << p->bonds.size() << "\n";
	for(int i1 = 1; i1 < N; i1++){
		for(int i2 = 0; i2 < i1; i2++){
			if(!p->bonds[i1]->outside || !p->bonds[i2]->outside){
				V += bondedForce(p, i1);
				V += bondedForce(p, i2);
				V += threeBodyForce(p, i1, i2);
			}
		}
	}
	return V;
}

void Potential::updateR(Particle *p, double dt){
	double dr;
	for(int d = 0; d < 3; d++){
		dr = dt*p->v[d];
		p->r[d] = p->r[d] + dr;
		p->rReal[d] = p->rReal[d] + dr;
	}
}

void Potential::updateV(Particle *p){
	double ff = fFactor[p->type];
	for(int d = 0; d < 3; d++){
		p->v[d] = p->v[d] + ff*p->F[d];
	}
}

void Potential::createForceFactor(double dt){
	for(int i = 0; i < types; i++){
		fFactor[i] = 0.5*dt/m[i];
	}
}

void Potential::createVelocityFactor(double dt, double zeta){
	for(int i = 0; i < types; i++){
		vFactor[i] = 1/(1 + m[i]*fFactor[i]*zeta);
	}
}

void Potential::updateVnh1(Particle *p, double zeta){
	double ff = fFactor[p->type];
	for(int d = 0; d < 3; d++){
		p->v[d] = p->v[d] + ff*(p->F[d] - zeta*p->v[d]*m[p->type]);
	}
}

void Potential::updateVnh2(Particle *p){
	double ff = fFactor[p->type];
	double vf = vFactor[p->type];
	for(int d = 0; d < 3; d++){
		p->v[d] = vf*(p->v[d] + ff*p->F[d]);
	}
}