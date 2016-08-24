#ifndef PARTICLE_H
#define PARTICLE_H
#include "IntList.h"
#include <string>
#include <vector>

using std::vector;

class ParticleContainer;

class Particle{
	public:
		double* r;
		double* v;
		double* F;
		double* rReal;
		int id;
		int type;
		bool outside;
		//number of Si and H particles neighbours
		double nSi;
		double nH;
	public:
		Particle(double* data, int id, int type);
		Particle(double* data, int id, int type, bool tmp);
		~Particle();
		void resetForces();
		void resetLists();
		void monitor();
		vector<Particle*> neighbours;
		vector<Particle*> bonds;
		bool monitored;
		vector<double> forces2B;
		vector<int> index2B;
		// vector<double> forces3B;
		// vector<int> index3B;
		double dot;
};

#endif