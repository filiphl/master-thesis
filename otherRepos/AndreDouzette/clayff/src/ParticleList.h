#ifndef PARTICLELIST_H
#define PARTICLELIST_H
#include "Particle.h"
#include "ParticleContainer.h"
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include <mpi.h>
#include <vector>

using std::vector;

class ParticleList{
	public:
		ParticleList(ParticleContainer* container);
		~ParticleList();
		void add(Particle* p);
		void remove(int id);
		Particle* getFirst();
		Particle* getNext();
		bool hasNext();
		int getLength();
		void makeArrays(int** stats, double** data, int* nStats, int* nData);
		void makeBoundaryArrays(int** stats, double** data, int* nStats, int* nData);
		void reset();
		void add(ParticleList* list);
		void deleteParticles();
		void save();
		void resume();
		void makeWeightArrays(int** stats, double** data, int* nStats, int* nData, int typeO);
		int length;
		int getLast();
		void makeArray();
		Particle** array;
	private:
		class Node{
			public:
				int val;
				Node* next;
				Node(int val);
				void add(int val);
				bool remove(int id);
				void deleteAll();
		};
	private:
		Node* iterator;
		Node* saved;
		Node* first;
		ParticleContainer* container;
	public:
		bool outside;
};

#endif