#include "Particle.h"
#include <iostream>
//Need to be included here to circumvent inclusion recursion in header file.
#include "ParticleContainer.h"
using namespace std;

Particle::Particle(double* data, int id, int type){
	r = &(data[0]);
	v = &(data[3]);
	F = &(data[6]);
	rReal = &(data[9]);
	nH = 0;
	nSi = 0;
	this->id = id;
	//Change this to a pointer to the array element in container?
	this->type = type;
	outside = false;
	//Reserving space to neighbour and bond lists. Increase necessary?
	neighbours.reserve(100);
	bonds.reserve(10);
	monitored = false;
}

Particle::Particle(double* data, int id, int type, bool tmp){
	//Used for initializing boundary particles
	r = &(data[0]);
	v = &(data[3]);
	F = &(data[6]);
	rReal = &(data[9]);
	nH = 0;
	nSi = 0;
	this->id = id;
	//Change this to a pointer to the array element in container?
	this->type = type;
	//Reserving space to bond lists. No neighbours required for boundary particles
	bonds.reserve(10);
	outside = true;
	monitored = false;
}

void Particle::monitor(){
	monitored = true;
	forces2B.reserve(300);
	index2B.reserve(100);
	// forces3B.reserve(30);
	// index3B.reserve(30);
}

Particle::~Particle(){
	resetLists();
}

void Particle::resetForces(){
	for(int d = 0; d < 3; d++){
		F[d] = 0;
	}
}

void Particle::resetLists(){
	neighbours.clear();
	bonds.clear();
	forces2B.clear();
	index2B.clear();
}