#include "ParticleList.h"

//TO DO: Implement proper iterator

ParticleList::ParticleList(ParticleContainer* container){
	this->container = container;
	iterator = 0;
	first = 0;
	length = 0;
	array = new Particle*[750];
}

ParticleList::~ParticleList(){
	if(first != 0){
		reset();
	}
	delete[] array;
}

void ParticleList::makeArray(){
	int i = 0;
	for(Node* n = first; n != 0; n = n->next){
		array[i] = container->particles[n->val];
		i++;
	}
}

void ParticleList::add(Particle* p){
	if(first == 0){
		first = new Node(p->id);
	}else{
		//Add new node to start of list
		Node* tmp = first;
		first = new Node(p->id);
		first->next = tmp;
	}
	length += 1;
}

void ParticleList::add(ParticleList* list){
	for(Particle* p = list->getFirst(); list->hasNext(); p = list->getNext()){
		add(p);
	}
}

void ParticleList::remove(int id){
	if(first == 0){
		return;
	}
	if(first->val == id){
		Node* tmp = first;
		first = first->next;
		delete tmp;
	}else{
		for(Node* n = first; n->next != 0; n = n->next){
			if(n->next->val == id){
				Node* tmp = n->next;
				n->next = n->next->next;
				delete tmp;
				break;
			}
		}
	}
	length += -1;
}

Particle* ParticleList::getFirst(){
	iterator = first;
	if(first != 0){
		return container->particles[first->val];
	}else{
		return 0;
	}
}

Particle* ParticleList::getNext(){
	iterator = iterator->next;
	if(iterator != 0){
		return container->particles[iterator->val];
	}else{
		return 0;
	}
}

void ParticleList::save(){
	saved = iterator;
}

void ParticleList::resume(){
	iterator = saved;
}

bool ParticleList::hasNext(){
	return iterator != 0;
}

int ParticleList::getLength(){
	return length;
}

void ParticleList::makeArrays(int** stats, double** data, int* nStats, int* nData){
	*nStats = length;
	*nData = 9*length;
	int* tmpStats = new int[*nStats];
	double* tmpData = new double[*nData];
	Particle* p;
	int ip = 0;
	for(Node* node = first; node != 0; node = node->next){
		p = container->get(node->val);
		//Add particle stats to array
		tmpStats[ip] = p->id;
		//Add particle data to arary
		tmpData[9*ip]     = p->r[0];
		tmpData[9*ip + 1] = p->r[1];
		tmpData[9*ip + 2] = p->r[2];
		tmpData[9*ip + 3] = p->v[0];
		tmpData[9*ip + 4] = p->v[1];
		tmpData[9*ip + 5] = p->v[2];
		tmpData[9*ip + 6] = p->rReal[0];
		tmpData[9*ip + 7] = p->rReal[1];
		tmpData[9*ip + 8] = p->rReal[2];
		ip++;
	}
	*data = tmpData;
	*stats = tmpStats;
}

//Same as above, without veleocities, bonds are comunicated as well(used for 3BF)
void ParticleList::makeBoundaryArrays(int** stats, double** data, int* nStats, int* nData){
	*nStats = length;
	*nData = 3*length;
	int* tmpStats = new int[*nStats];
	double* tmpData = new double[*nData];
	Particle* p;
	int ip = 0;
	for(Node* node = first; node != 0; node = node->next){
		p = container->get(node->val);
		//Add particle stats to array
		tmpStats[ip] = p->id;
		//Add particle data to arary
		tmpData[3*ip    ] = p->r[0];
		tmpData[3*ip + 1] = p->r[1];
		tmpData[3*ip + 2] = p->r[2];
		ip++;
	}
	*data = tmpData;
	*stats = tmpStats;
}

//make more efficient!
void ParticleList::makeWeightArrays(int** stats, double** data, int* nStats, int* nData, int typeO){
	Particle* p;
	int oxygens = 0;
	for(Node* n = first; n != 0; n = n->next){
		p = container->get(n->val);
		if(p->type == typeO){
			oxygens += 1;
		}
	}
	*nStats = oxygens;
	*nData = 2*oxygens;
	int* tmpStats = new int[*nStats];
	double* tmpData = new double[*nData];
	int i = 0;
	for(Node* n = first; n != 0; n = n->next){
		p = container->get(n->val);
		if(p->type == typeO){
			tmpStats[i] = p->id;
			tmpData[2*i] = p->nSi;
			tmpData[2*i + 1] = p->nH;
			i += 1;
		}
	}
	*data = tmpData;
	*stats = tmpStats;
}

void ParticleList::deleteParticles(){
	//Delete nodes and particles
	Node* nn;
	for(Node* n = first; n != 0; n = nn){
		nn = n->next;
		container->del(n->val);
		delete n;
	}
	first = 0;
	length = 0;
}

void ParticleList::reset(){
	//Delete nodes
	Node* nn;
	for(Node* n = first; n != 0; n = nn){
		nn = n->next;
		delete n;
	}
	first = 0;
	length = 0;
}

////////////////
//NODE METHODS//
////////////////

void ParticleList::Node::deleteAll(){
	if(next != 0){
		next->deleteAll();
	}
	delete this;
}

ParticleList::Node::Node(int val){
	this->val = val;
	next = 0;
}

void ParticleList::Node::add(int val){
	if(next == 0){
		next = new Node(val);
	}else{
		next->add(val);
	}
}

bool ParticleList::Node::remove(int id){
	if(next == 0){
		return false;
	}
	if(next->val == id){
		Node* tmp = next;
		next = next->next;
		delete tmp;
		return true;
	}
	return next->remove(id);
}

// #include "ParticleList.h"

// ParticleList::ParticleList(ParticleContainer* container){
// 	this->container = container;
// 	iterator = 0;
// 	first = 0;
// 	length = 0;
// 	array = new Particle*[250];
// }

// ParticleList::~ParticleList(){
// 	if(first != 0){
// 		reset();
// 	}
// 	delete[] array;
// }

// void ParticleList::makeArray(){
// 	int i = 0;
// 	for(Node* n = first; n != 0; n = n->next){
// 		array[i] = container->particles[n->val];
// 		i++;
// 	}
// }

// void ParticleList::add(Particle* p){
// 	if(first == 0){
// 		first = new Node(p->id);
// 	}else{
// 		//Add new node to start of list
// 		Node* tmp = first;
// 		first = new Node(p->id);
// 		first->next = tmp;
// 	}
// 	length += 1;
// }

// void ParticleList::add(ParticleList* list){
// 	for(Particle* p = list->getFirst(); list->hasNext(); p = list->getNext()){
// 		add(p);
// 	}
// }

// void ParticleList::remove(int id){
// 	if(first == 0){
// 		return;
// 	}
// 	if(first->val == id){
// 		Node* tmp = first;
// 		first = first->next;
// 		delete tmp;
// 	}else{
// 		for(Node* n = first; n->next != 0; n = n->next){
// 			if(n->next->val == id){
// 				Node* tmp = n->next;
// 				n->next = n->next->next;
// 				delete tmp;
// 				break;
// 			}
// 		}
// 	}
// 	length += -1;
// }

// Particle* ParticleList::getFirst(){
// 	iterator = first;
// 	if(first != 0){
// 		return container->particles[first->val];
// 	}else{
// 		return 0;
// 	}
// }

// Particle* ParticleList::getNext(){
// 	iterator = iterator->next;
// 	if(iterator != 0){
// 		return container->particles[iterator->val];
// 	}else{
// 		return 0;
// 	}
// }

// void ParticleList::save(){
// 	saved = iterator;
// }

// void ParticleList::resume(){
// 	iterator = saved;
// }

// bool ParticleList::hasNext(){
// 	return iterator != 0;
// }

// int ParticleList::getLength(){
// 	return length;
// }

// void ParticleList::makeArrays(int** stats, double** data, int* nStats, int* nData){
// 	*nStats = length;
// 	*nData = 9*length;
// 	int* tmpStats = new int[*nStats];
// 	double* tmpData = new double[*nData];
// 	Particle* p;
// 	int ip = 0;
// 	for(Node* node = first; node != 0; node = node->next){
// 		p = container->get(node->val);
// 		//Add particle stats to array
// 		tmpStats[ip] = p->id;
// 		//Add particle data to arary
// 		tmpData[9*ip]     = p->r[0];
// 		tmpData[9*ip + 1] = p->r[1];
// 		tmpData[9*ip + 2] = p->r[2];
// 		tmpData[9*ip + 3] = p->v[0];
// 		tmpData[9*ip + 4] = p->v[1];
// 		tmpData[9*ip + 5] = p->v[2];
// 		tmpData[9*ip + 6] = p->rReal[0];
// 		tmpData[9*ip + 7] = p->rReal[1];
// 		tmpData[9*ip + 8] = p->rReal[2];
// 		ip++;
// 	}
// 	*data = tmpData;
// 	*stats = tmpStats;
// }

// //Same as above, without veleocities, bonds are comunicated as well(used for 3BF)
// void ParticleList::makeBoundaryArrays(int** stats, double** data, int* nStats, int* nData){
// 	*nStats = length;
// 	*nData = 3*length;
// 	int* tmpStats = new int[*nStats];
// 	double* tmpData = new double[*nData];
// 	Particle* p;
// 	int ip = 0;
// 	for(Node* node = first; node != 0; node = node->next){
// 		p = container->get(node->val);
// 		//Add particle stats to array
// 		tmpStats[ip] = p->id;
// 		//Add particle data to arary
// 		tmpData[3*ip    ] = p->r[0];
// 		tmpData[3*ip + 1] = p->r[1];
// 		tmpData[3*ip + 2] = p->r[2];
// 		ip++;
// 	}
// 	*data = tmpData;
// 	*stats = tmpStats;
// }

// //make more efficient!
// void ParticleList::makeWeightArrays(int** stats, double** data, int* nStats, int* nData, int typeO){
// 	Particle* p;
// 	int oxygens = 0;
// 	for(Node* n = first; n != 0; n = n->next){
// 		p = container->get(n->val);
// 		if(p->type == typeO){
// 			oxygens += 1;
// 		}
// 	}
// 	*nStats = oxygens;
// 	*nData = 2*oxygens;
// 	int* tmpStats = new int[*nStats];
// 	double* tmpData = new double[*nData];
// 	int i = 0;
// 	for(Node* n = first; n != 0; n = n->next){
// 		p = container->get(n->val);
// 		if(p->type == typeO){
// 			tmpStats[i] = p->id;
// 			tmpData[2*i] = p->nSi;
// 			tmpData[2*i + 1] = p->nH;
// 			i += 1;
// 		}
// 	}
// 	*data = tmpData;
// 	*stats = tmpStats;
// }

// void ParticleList::deleteParticles(){
// 	//Delete nodes and particles
// 	Node* nn;
// 	for(Node* n = first; n != 0; n = nn){
// 		nn = n->next;
// 		container->del(n->val);
// 		delete n;
// 	}
// 	first = 0;
// 	length = 0;
// }

// void ParticleList::reset(){
// 	//Delete nodes
// 	Node* nn;
// 	for(Node* n = first; n != 0; n = nn){
// 		nn = n->next;
// 		delete n;
// 	}
// 	first = 0;
// 	length = 0;
// }

// ////////////////
// //NODE METHODS//
// ////////////////

// void ParticleList::Node::deleteAll(){
// 	if(next != 0){
// 		next->deleteAll();
// 	}
// 	delete this;
// }

// ParticleList::Node::Node(int val){
// 	this->val = val;
// 	next = 0;
// }

// void ParticleList::Node::add(int val){
// 	if(next == 0){
// 		next = new Node(val);
// 	}else{
// 		next->add(val);
// 	}
// }

// bool ParticleList::Node::remove(int id){
// 	if(next == 0){
// 		return false;
// 	}
// 	if(next->val == id){
// 		Node* tmp = next;
// 		next = next->next;
// 		delete tmp;
// 		return true;
// 	}
// 	return next->remove(id);
// }