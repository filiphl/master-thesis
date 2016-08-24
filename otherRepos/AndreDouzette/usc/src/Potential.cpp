#include "Potential.h"
#include <iostream>
#include <mpi.h>
using namespace std;

//The USC model
Potential::Potential(){
	tau = 6.2831853071795864769252867665590057683943387987502116; //truncate
	taui = 1/tau;
	pi = 0.5*tau;
	//Correction factor for wrong units in paper
	double bohr = 0.52917721092;
	double atomicEnergyUnit = 27.21138;
	double correctionFactor = bohr*atomicEnergyUnit;
	
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
	//Charge 1e: THIS IS NOT NEEDED
	q0 = 1;
	
	//Indexes types
	typeSi = 0;
	typeO = 1;
	typeH = 2;
	typeOs = 1;
	typeOh = 3;
	
	//Number of particle types
	types = 4;
	
	//2BF cutoff
	rc  = 5.5;
	
	correctionFactor = correctionFactor/E0;
	
	//Particle properties
	double* z = new double[types];
	m = new double[types];
	fFactor = new double[types];
	vFactor = new double[types];
	//Interpolation parameters
	Rb = new double[types];
	Db = new double[types];
	Dbi = new double[types];
	rInner2 = new double[types];
	rOuter2 = new double[types];
	//Two body interactions
	Z = new double*[types];
	H  = new double*[types];
	D  = new double*[types];
	W  = new double*[types];
	eta = new double*[types];
	//Shielding length
	double **r1s = new double*[types];
	double **r4s = new double*[types];
	r1si = new double*[types];
	r4si = new double*[types];
	for(int i = 0; i < types; i++){
		Z[i] = new double[types];
		H[i]  = new double[types];
		D[i]  = new double[types];
		W[i]  = new double[types];
		eta[i] = new double[types];
		r1s[i] = new double[types];
		r4s[i] = new double[types];
		r1si[i] = new double[types];
		r4si[i] = new double[types];
		//This is probably not needed
		for(int j = 0; j < types; j++){
			Z[i][j] = 0;
			H[i][j] = 0;
			D[i][j] = 0;
			W[i][j] = 0;
			eta[i][j] = 0;
			r1s[i][j] = 0;
			r4s[i][j] = 0;
			r1si[i][j] = 0;
			r4si[i][j] = 0;
		}
	}
	
	//Three body interactions
	double*** theta0 = new double**[types];
	B  = new double**[types];
	c0 = new double**[types];
	xi = new double*[types];
	r0 = new double*[types];
	r02 = new double*[types];
	
	for(int i = 0; i < types; i++){
		B[i]  = new double*[types];
		theta0[i] = new double*[types];
		c0[i] = new double*[types];
		xi[i] = new double[types];
		r0[i] = new double[types];
		r02[i] = new double[types];
		for(int j = 0; j < types; j++){
			B[i][j]  = new double[types];
			theta0[i][j] = new double[types];
			c0[i][j] = new double[types];
			r0[i][j] = 0;
			xi[i][j] = 0;
			for(int k = 0; k < types; k++){
				B[i][j][k] = 0;
				theta0[i][j][k] = 0;
				c0[i][j][k] = 0;
			}
		}
	}
	
	
	//Screening of coulomb
	r1s[typeSi][typeSi] = 4.43;
	r1s[typeSi][typeO]  = 4.43;
	r1s[typeOs][typeOs] = 4.43;
	r1s[typeO][typeH]   = 4.43;
	r1s[typeH][typeH]   = 4.43;
	r1s[typeOh][typeOh] = 4.43;
	r1s[typeSi][typeH]  = r1s[typeH][typeH];
	//Screening of dipole
	r4s[typeSi][typeSi] = 2.50;
	r4s[typeSi][typeO]  = 2.50;
	r4s[typeOs][typeOs] = 2.50;
	r4s[typeO][typeH]   = 1.51113;
	r4s[typeH][typeH]   = 2.50;
	r4s[typeOh][typeOh] = 2.50;

	//Charge
	z[typeSi] = 1.2;
	z[typeOs] = -0.6;
	z[typeH]  = 0.32983;
	z[typeOh] = -0.65966;
	//Mass
	m[typeSi] = 28.0855;
	m[typeO]  = 15.9994;
	m[typeH]  = 1.00794;
	m[typeOh] = 1; //This is never used
	//Steric repulsion
	eta[typeSi][typeSi] = 11;
	eta[typeSi][typeO]  = 9;
	eta[typeOs][typeOs] = 7;
	eta[typeH][typeH]   = 9;
	eta[typeO][typeH]   = 9;
	eta[typeOh][typeOh] = 9;
	H[typeSi][typeSi] = 0.39246;
	H[typeSi][typeO]  = 78.3143;
	H[typeOs][typeOs] = 355.5263;
	H[typeO][typeH]   = 0.61437;
	H[typeOh][typeOh] = 1965.88;
	//Charge-dipole Strength
	D[typeSi][typeO]  = 3.456;
	D[typeOs][typeOs] = 1.728;
	D[typeO][typeH]   = 0.2611;
	D[typeOh][typeOh] = 2.0887;
	//Van der waals interaction
	W[typeOh][typeOh] = 10.0;
	
	//3BF Strength
	B[typeO][typeSi][typeO]  = 4.993;
	B[typeSi][typeO][typeSi] = 19.972;
	B[typeH][typeO][typeH]   = 52.9333;
	B[typeSi][typeO][typeH]  = 36.45265;
	//Equilibrium angle
	theta0[typeO][typeSi][typeO]  = 109.47;
	theta0[typeSi][typeO][typeSi] = 141.0;
	theta0[typeH][typeO][typeH]   = 97.9476;
	theta0[typeSi][typeO][typeH]  = 110.0;
	//Screening factor
	xi[typeSi][typeO] = 1.0;
	xi[typeO][typeH]   = 0.75;
	//3BF cutoff
	r0[typeSi][typeO] = 2.6;
	r0[typeO][typeH]  = 1.4;
	
	//Oxygen number interpolation parameters
	Rb[typeSi] = 2.0;
	Rb[typeH]  = 1.4;
	Db[typeSi] = 0.3;
	Db[typeH]  = 0.3;
	
	
	//Symmetric
	for(int i = 0; i < types; i++){
		for(int j = 0; j < i; j++){
			eta[i][j] = eta[j][i];
			H[i][j] = H[j][i];
			D[i][j] = D[j][i];
			W[i][j] = W[j][i];
			r1s[i][j] = r1s[j][i];
			r4s[i][j] = r4s[j][i];
			r0[i][j] = r0[j][i];
			xi[i][j] = xi[j][i];
			for(int k = 0; k < types; k++){
				B[i][k][j] = B[j][k][i];
				theta0[i][k][j] = theta0[j][k][i];
			}
		}
	}
	
	
	//Derived values
	for(int i = 0; i < types; i++){
		Dbi[i] = 1/Db[i];
		for(int j = 0; j < types; j++){
			r1si[i][j] = 1/r1s[i][j];
			r4si[i][j] = 1/r4s[i][j];
			r02[i][j] = r0[i][j]*r0[i][j];
			for(int k = 0; k < types; k++){
				c0[i][j][k] = cos(theta0[i][j][k]*pi/180);
			}
			delete[] theta0[i][j];
		}
		delete[] r1s[i];
		delete[] r4s[i];
		delete[] theta0[i];
	}
	delete[] r1s;
	delete[] r4s;
	delete[] theta0;
	
	
	double L0eta;
	double L02 = L0*L0;
	double L04 = L02*L02;
	double L06 = L04*L02;
	double rInner, rOuter;
	
	//Normalize to numerical values
	for(int i = 0; i < types; i++){
		z[i] = z[i]/q0;
		m[i] = m[i]/m0;
		Rb[i] = Rb[i]/L0;
		Db[i] = Db[i]/L0;
		rInner = Rb[i] + Db[i];
		rOuter = Rb[i] - Db[i];
		rInner2[i] = rInner*rInner;
		rOuter2[i] = rOuter*rOuter;
		for(int j = 0; j < types; j++){
			//calculate L0^eta
			L0eta = L0;
			for(int k = 0; k < eta[i][j] - 1; k++){
				L0eta = L0eta*L0;
			}
			Z[i][j] = z[i]*z[j]*correctionFactor;
			H[i][j] = H[i][j]/(L0eta*E0);
			D[i][j] = D[i][j]*correctionFactor/L04;
			W[i][j] = W[i][j]/(L06*E0);
			r1si[i][j] = r1si[i][j]*L0;
			r4si[i][j] = r4si[i][j]*L0;
			r0[i][j] = r0[i][j]/L0;
			r02[i][j] = r0[i][j]*r0[i][j];
			xi[i][j] = xi[i][j]/L0;
			for(int k = 0; k < types; k++){
				B[i][j][k] = B[i][j][k]/E0;
			}
		}
	}
	rc = rc/L0;
	rc2 = rc*rc;
	delete[] z;
	
	//Some fixes for O-H interaction
	Z[typeO][typeH] = z[typeOh]*z[typeH]*correctionFactor;
	
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
		for(int j = 0; j <= i; j++){
			Tabulate2BPotential(i, j);
		}
	}
	//Same array for the symetric part
	for(int i = 0; i < types; i++){
		for(int j = i + 1; j < types; j++){
			fTable[i][j] = fTable[j][i];
			VTable[i][j] = VTable[j][i];
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

//Remove? used for testing purposes
void Potential::writePotentialFile(string name, double* r2,  double* f, double* V, int n){
	int pres = 10;
	ofstream potentialFile;
	potentialFile.open(name.c_str());
	if(potentialFile.is_open()){
		cout << "writing " << name << "\n";
		for(int i = 0; i < n; i++){
			potentialFile << setprecision(pres) << r2[i] << " ";
			potentialFile << setprecision(pres) << f[i] << " ";
			potentialFile << setprecision(pres) << V[i] << "\n";
		}
		potentialFile.close();
	}else{
		cout << "Error writing to file " << name << "\n";
	}
}

void Potential::print(){
	for(int i = 0; i < types; i++){
		for(int j = 0; j < types; j++){
			cout << "\n\n\n\n" << i << ", " << j << "\n";
			cout << "P " << H[i][j] << "\t\t" << Z[i][j] << "\t\t" << D[i][j] << "\t\t" << W[i][j] << "\n";
			cout << "P " << r1si[i][j] << "\t\t" << r4si[i][j] << "\t\t" << eta[i][j] << "\n";
			for(int k = 0; k < tableSize; k++){
				cout << r2Table[k] << "\t\t" << fTable[i][j][k] << "\t\t" << VTable[i][j][k] << "\n";
			}
		}
	}
}

void Potential::Tabulate2BPotential(int type1, int type2){
	double r, r2, r4, r6, r8;
	double reta, tmp;
	double* rv = new double[3];
	double f1, f2, f3, f4;
	double V1, V2, V3, V4;
	double Vc, fc;
	int n = tableSize;
	//Importing correct particle parameters
	double H = this->H[type1][type2];
	double D = this->D[type1][type2];
	double W = this->W[type1][type2];
	double eta = this->eta[type1][type2];
	double r1si = this->r1si[type1][type2];
	double r4si = this->r4si[type1][type2];
	double Z = this->Z[type1][type2];
	//Calculating potential in nodes given in rTable
	fTable[type1][type2] = new double[n];
	VTable[type1][type2] = new double[n];
	double* f = fTable[type1][type2];
	double* V = VTable[type1][type2];
	
	for(int i = 0; i < n; i++){
		r2 = r2Table[i];
		r = sqrt(r2);
		r4 = r2*r2;
		r6 = r4*r2;
		r8 = r4*r4;
		//Calculating reta
		reta = r;
		for(int j = 0; j < eta - 1; j++){
			reta = reta*r;
		}
		//Force calculation
		if(H == 0){
			V1 = 0;
			f1 = 0;
		}else{
			V1 = H/reta;
			f1 = V1*eta/r2;
		}
		if(Z == 0){
			V2 = 0;
			f2 = 0;
		}else{
			tmp = r*r1si;
			V2 = Z/r*exp(-tmp);
			f2 = V2*(tmp + 1)/r2;
		}
		if(D == 0){
			V3 = 0;
			f3 = 0;
		}else{
			tmp = r*r4si;
			V3 = -D/(2*r4)*exp(-tmp);
			f3 = V3*(tmp + 4)/r2;
		}
		if(W == 0){
			V4 = 0;
			f4 = 0;
		}else{
			V4 = -W/r6;
			f4 = V4*6/r2;
		}
		//Adding up force and potential components
		f[i] = f1 + f2 + f3 + f4;
		V[i] = V1 + V2 + V3 + V4;
	}
	//Truncating
	Vc = V[n - 1];
	fc = f[n - 1];
	for(int i = 0; i < n; i++){
		r = sqrt(r2Table[i]);
		V[i] = V[i] - Vc + (r - rc)*rc*fc; //f is stored as F/r
		f[i] = f[i] - fc*rc/r;
	}
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
	//Covnerting to vectors and adding forces to particles(f = |F|/r)
	for(int d = 0; d < 3; d++){
		Fv = f*rv[d];
		p->F[d] += -Fv;
		q->F[d] +=  Fv;
	}
	
	// if(p->outside || q->outside){
	if(q->outside){
		return 0.5*V;
	}
	return V;
}

double Potential::twoBodyOxygenForce(Particle *p, int bond){
	Particle* q = p->neighbours[bond];
	double t;
	int n;
	double *ft, *Vt;
	double fs, fh, Vs, Vh, f, V, Fv;
	int nSi;
	int nH;
	//Find interval in which to interpolate
	double r2 = 0;
	for(int d = 0; d < 3; d++){
		rv[d] = q->r[d] - p->r[d];
		r2 += rv[d]*rv[d];
	}
	
	r2 = tableIncrementInv*r2;
	n = (int)r2;
	//Interpolating scheme for table lookup
	t = r2 - n;
	ft = fTable[typeOs][typeOs];
	Vt = VTable[typeOs][typeOs];
	fh = ft[n]*(1 - t) + t*ft[n + 1];
	Vh = Vt[n]*(1 - t) + t*Vt[n + 1];
	ft = fTable[typeOh][typeOh];
	Vt = VTable[typeOh][typeOh];
	fs = ft[n]*(1 - t) + t*ft[n + 1];
	Vs = Vt[n]*(1 - t) + t*Vt[n + 1];
	//Interpolation between Oh and Os
	nSi = p->nSi + q->nSi;
	nH = p->nH + q->nH;
	t = nH/((double)(nH + nSi));
	V = Vs*(1 - t) + Vh*t;
	f = fs*(1 - t) + fh*t;
	//Covnerting to vectors and adding forces to particles(f = |F|/r)
	
	for(int d = 0; d < 3; d++){
		Fv = f*rv[d];
		p->F[d] += -Fv;
		q->F[d] +=  Fv;
	}
	
	// if(p->outside || q->outside){
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

double Potential::twoBodyOxygenForceMonitoring(Particle *p, int bond){
	Particle* q = p->neighbours[bond];
	double t;
	int n;
	double *ft, *Vt;
	double fs, fh, Vs, Vh, f, V, Fv;
	int nSi;
	int nH;
	//Find interval in which to interpolate
	double r2 = 0;
	for(int d = 0; d < 3; d++){
		rv[d] = q->r[d] - p->r[d];
		r2 += rv[d]*rv[d];
	}
	
	r2 = tableIncrementInv*r2;
	n = (int)r2;
	//Interpolating scheme for table lookup
	t = r2 - n;
	ft = fTable[typeOs][typeOs];
	Vt = VTable[typeOs][typeOs];
	fh = ft[n]*(1 - t) + t*ft[n + 1];
	Vh = Vt[n]*(1 - t) + t*Vt[n + 1];
	ft = fTable[typeOh][typeOh];
	Vt = VTable[typeOh][typeOh];
	fs = ft[n]*(1 - t) + t*ft[n + 1];
	Vs = Vt[n]*(1 - t) + t*Vt[n + 1];
	//Interpolation between Oh and Os
	nSi = p->nSi + q->nSi;
	nH = p->nH + q->nH;
	t = nH/((double)(nH + nSi));
	V = Vs*(1 - t) + Vh*t;
	f = fs*(1 - t) + fh*t;
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
	// return 0;
	Particle* q1 = p->bonds[bond1];
	Particle* q2 = p->bonds[bond2];
	double B = this->B[q1->type][p->type][q2->type];
	if(B == 0) return 0;
	
	//Retrieving parameters
	double c0 = this->c0[q1->type][p->type][q2->type];
	double xi1 = this->xi[p->type][q1->type];
	double xi2 = this->xi[p->type][q2->type];
	double r01 = this->r0[p->type][q1->type];
	double r02 = this->r0[p->type][q2->type];
	double dc;
	double r1i, r2i, dr1i, dr2i;
	double Vt, a1, b1, b2;
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
	dr1i = 1/(r1 - r01);
	dr2i = 1/(r2 - r02);
	r1i = 1/r1;
	r2i = 1/r2;
	c = c*r1i*r2i;
	dc = c - c0;
	Vt = B*dc*exp(xi1*dr1i + xi2*dr2i);
	//Force matrix components
	a1 = Vt*r1i*(2*c*r1i + dc*xi1*dr1i*dr1i);
	b1 = -2*Vt*r1i*r2i;
	b2 = Vt*r2i*(2*c*r2i + dc*xi2*dr2i*dr2i);
	
	//Add collection of monitored particle forces?
	for(int d = 0; d < 3; d++){
		F1 = a1*rv1[d] + b1*rv2[d];
		F2 = b1*rv1[d] + b2*rv2[d];
		q1->F[d] += F1;
		q2->F[d] += F2;
		p->F[d] += -(F1 + F2);
	}
	//Potential energy calculation
	return Vt*dc*(1 - oneThird*(p->outside + q1->outside + q2->outside));
}

double Potential::force2(Particle* p){
	double V = 0;
	int N = p->neighbours.size();
	if(p->type == typeO){
		int nSi, nH;
		//Have a variable instead of going through p->neighbours???
		for(int i = 0; i < N; i++){
			if(p->neighbours[i]->type == typeO){
				nH = p->nH + p->neighbours[i]->nH;
				nSi = p->nSi + p->neighbours[i]->nSi;
				if(nH == 0){
					V += twoBodyForce(p, i);
				}else if(nSi == 0){
					p->type = typeOh;
					p->neighbours[i]->type = typeOh;
					V += twoBodyForce(p, i);
					p->type = typeO;
					p->neighbours[i]->type = typeO;
				}else{
					V += twoBodyOxygenForce(p, i);
				}
			}else{
				V += twoBodyForce(p, i);
			}
		}
	}else{
		for(int i = 0; i < N; i++){
			V += twoBodyForce(p, i);
		}
	}
	return V;
}

double Potential::force2Monitoring(Particle* p){
	double V = 0;
	int N = p->neighbours.size();
	if(p->type == typeO){
		int nSi, nH;
		//Have a variable instead of going through p->neighbours???
		for(int i = 0; i < N; i++){
			if(p->neighbours[i]->type == typeO){
				nH = p->nH + p->neighbours[i]->nH;
				nSi = p->nSi + p->neighbours[i]->nSi;
				if(nH == 0){
					V += twoBodyForceMonitoring(p, i);
				}else if(nSi == 0){
					p->type = typeOh;
					p->neighbours[i]->type = typeOh;
					V += twoBodyForceMonitoring(p, i);
					p->type = typeO;
					p->neighbours[i]->type = typeO;
				}else{
				V += twoBodyOxygenForceMonitoring(p, i);
				}
			}else{
				V += twoBodyForceMonitoring(p, i);
			}
		}
		// if(p->monitored){
		// 	cout << p->index2B.size() << "\n";
		// }
	}else{
		for(int i = 0; i < N; i++){
			V += twoBodyForceMonitoring(p, i);
		}
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
	for(int i1 = 1; i1 < N; i1++){
		for(int i2 = 0; i2 < i1; i2++){
			if(!p->bonds[i1]->outside || !p->bonds[i2]->outside){
				V += threeBodyForce(p, i1, i2);
			}
		}
	}
	return V;
}

double Potential::numberInterpolate(Particle *p, Particle *pb, double r2){
	if(r2 >= rOuter2[pb->type]){
		return 0;
	}else if(r2 < rInner2[pb->type]){
		return 1;
	}else{
		double A = Dbi[pb->type]*(sqrt(r2) - Rb[pb->type]) + 1;
		return 1 - 0.5*A + taui*sin(pi*A);
	}
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
		vFactor[i] = 1/(1 + 0.5*dt*zeta);
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