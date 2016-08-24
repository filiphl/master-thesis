#include "ParticleContainer.h"
#include <iostream>
using namespace std;
//For easy access to particles, when we know the id
//TO DO: Add array for position, velocity and force storage

// Indexed array for particle pointer storrage
// Store positions, velocities and forces as an array-list?

ParticleContainer::ParticleContainer(int N){
	this->N = N;
	particles = new Particle*[N];
	mark = new bool[N];
	data = new double[12*N];
	types = new int[N];
	monitored = new bool[N];
	for(int i = 0; i < N; i++){
		particles[i] = 0;
		types[i] = 0;
		mark[i] = false;
		monitored[i] = false;
	}
	bonds = new int*[N];
	for(int i = 0; i < N; i++){
		bonds[i] = new int[2];
		bonds[i][0] = -1;
		bonds[i][1] = -1;
	}
}

void ParticleContainer::addBond(int source, int bond){
	//Hardcoded as only two bonds permitted
	if(bonds[source][0] == -1){
		bonds[source][0] = bond;
	}else{
		bonds[source][1] = bond;
	}
}

Particle* ParticleContainer::createParticle(double x, double y, double z, double vx, double vy, double vz, int id){
	int type = types[id];
	int i = 12*id;
	data[i] = x;
	data[i + 1] = y;
	data[i + 2] = z;
	data[i + 3] = vx;
	data[i + 4] = vy;
	data[i + 5] = vz;
	//Force
	data[i + 6] = 0;
	data[i + 7] = 0;
	data[i + 8] = 0;
	particles[id] = new Particle(&(data[i]), id, type);
	particles[id]->monitored = monitored[id];
	return particles[id];
}

Particle* ParticleContainer::createParticle(double x, double y, double z, int id){
	int type = types[id];
	int i = 12*id;
	data[i] = x;
	data[i + 1] = y;
	data[i + 2] = z;
	particles[id] = new Particle(&(data[i]), id, type, true);
	return particles[id];
}


int ParticleContainer::getNumberOfParticles(){
	int np = 0;
	for(int i = 0; i < N; i++){
		if(particles[i] != 0){
			np += 1;
		}
	}
	return np;
}

void ParticleContainer::add(Particle* p){
	particles[p->id] = p;
}

void ParticleContainer::remove(int id){
	particles[id] = 0;
}

void ParticleContainer::del(int id){
	if(particles[id] != 0){
		delete particles[id];
		particles[id] = 0;
	}
}

void ParticleContainer::delAll(){
	for(int i = 0; i < N; i++){
		if(particles[i] != 0){
			del(i);
		}
	}
}

Particle* ParticleContainer::get(int id){
	return particles[id];
}