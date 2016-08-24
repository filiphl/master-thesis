#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

#include "Particle.h"
#include <vector>
using std::vector;

class ParticleContainer{
	public:
		Particle** particles;
		bool* mark;
		bool* monitored;
		int N;
		double* data;
		int* types;
		int** bonds;
	public:
		ParticleContainer(int N);
		Particle* createParticle(double x, double y, double z, double vx, double vy, double vz, int id);
		Particle* createParticle(double x, double y, double z, int id);
		void add(Particle* p);
		void remove(int id);
		Particle* get(int id);
		void del(int id);
		void delAll();
		int getNumberOfParticles();
		vector<double> neighData;
		vector<Particle*> neighbours;
		void addBond(int, int);
};

#endif