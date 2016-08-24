#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

class NeighborList{
	public:
		NeigborLits();
		~NeighborList();
		void add(int id, double r2);
		void getFirst(int *idget, double *r2get, double *rget);
		void getNext(int *idget, double *r2get, double *rget);
		int getLength();
		void remove(int id);
		bool exists(int id);
		void save();
		void resume();
		void reset();
		
	private: class Node{
			private:
				ParticleInfo info;
			public:
				Node* next;
				Node(int id, double r2, double *r);
				void deleteAll();
				void add(int id, double r2, double *r);
				int remove(int id);
				void get(int *idget, double *r2get, double *rget);
				int get();
		};
	private:
		Node *iterator;
		Node *first;
		Node *last;
		int length;
		Node *saved;
};

#endif