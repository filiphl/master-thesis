#include "usc.h"

//Classes needed: particle, lists


//Initialization of parameters

void usc::init(){
	
}

//1 body force function

void usc::body1Force(){
	return;
}

//2 body force function
void usc::body2Force(Particle* p1, Particle* p2){
	double r2 = 0;
	double r, r4, r8;
	double rnu = 1;
	double* rv = new double[3];
	double F, F1, F2, F3, F4, tmpF;
	//Importing correct particle parameters
	double H = this->H[p1->type][p2->type];
	double D = this->D[p1->type][p2->type];
	double W = this->W[p1->type][p2->type];
	double nu = this->nu[p1->type][p2->type];
	double r1s = this->r1s[p1->type][p2->type];
	double r4s = this->r4s[p1->type][p2->type];
	double Z = this->Z[p1->type]*this->Z[p2->type];
	//Distance between particles
	for(int d = 0; d < 3; d++){
		rv[d] = p2->r[d] - p1->r[d];
		r2 += rv[d]*rv[d];
	}
	r = sqrt(r2);
	r4 = r2*r2;
	r8 = r4*r4;
	//Calculating rnu [find a better solution?]
	for(int n = 0; n < nu; n++){
		rnu = rnu*r;
	}
	//Force calculation
	F1 = nu*H/(rnu*r);
	F2 = Z*exp(-r/r1s)*(1/r1s + 1/r)/r2;
	F3 = -D/(2*r4*r)*exp(-r/r4s)*(1/r4s + 4/r);
	F4 = -6*W/r8;
	F = F1 + F2 + F3 + F4;
	//Applying force to particles
	for(d = 0; d < 3; d++){
		tmpF = F*rv[d];
		p1->F[d] +=  tmpF;
		p2->F[d] += -tmpF
	}
	delete[] rv;
}

//3 body force function
void usc::body3Force(Particle* p0, Particle* p1, Particle* p2){
	double r21 = 0;
	double r22 = 0;
	double c = 0;
	double r1, r2;
	double f1, f2;
	double ir1, ir2;
	double* rv1 = new double[3];
	double* rv2 = new double[3];
	for(int d = 0; d < 3; d++){
		rv1[d] = p1->r[d] - p0->r[d];
		rv2[d] = p2->r[d] - p0->r[d];
		r21 += rv1[d]*rv1[d];
		r22 += rv2[d]*rv2[d];
		c  += rv1[d]*rv2[d];
	}
	//Both particles must be within range
	if(r21 <= r0 || r22 <= r0){
		return;
	}
	r1 = sqrt(r21);
	r2 = sqrt(r22);
	c = c/(r1*r2);
	ir1 = 1/(r1 - r0);
	ir2 = 1/(r2 - r0);
	
	double dc = (c - c0);
	double Ve = B*exp(xi*(ir1 + ir2))
	double Fe = Ve*xi*dc*dc;
	double F1 = Fe*ir1*ir1;
	double F2 = Fe*ir2*ir2;
	
	ir1 = 1/r1;
	ir2 = 1/r2;
	for(int d = 0; d < 3; d++){
		f1 = F1*rv1[d]*ir1 - Ve*dc*ir1*(rv1[d]*ir1*c - rv2[d]*ir2);
		f2 = F2*rv2[d]*ir2 - Ve*dc*ir2*(rv2[d]*ir2*c - rv1[d]*ir1);
		p1->F[d] += f1;
		p2->F[d] += f2;
		p0->F[d] += -(f1 + f2);
	}
	
	delete[] rv1;
	delete[] rv2;
}

//Find bonds (if r < r0 then apply bond, else remove bond?)

//1. Init parameters


//2. Force calculation
