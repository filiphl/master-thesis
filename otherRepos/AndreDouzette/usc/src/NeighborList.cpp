#include "NeighborList.h"

using namespace std;

NeighborList::NeighborList(){
	iterator = 0;
	first = 0;
	last = 0;
	length = 0;
}

NeighborList::~NeighborList(){
	Node* nn;
	for(Node* n = first; n != 0; n = nn){
		nn = n->next;
		delete n;
	}
}

void NeighborList::add(int id, double r2, double *r){
	if(first == 0){
		first = new Node(id, r2, r);
		last = first;
	}else{
		last->add(id, r2, r);
		last = last ->next;
	}
	length += 1;
}

int NeighborList::getFirst(int *idget, double *r2get, double *rget){
	iterator = first;
	if(first != 0){
		first->get(idget, r2get, rget);
	}
}

bool NeighborList::hasNext(){
	return iterator != 0;
}

int NeighborList::getNext(int *idget, double *r2get, double *rget){
	iterator = iterator->next;
	if(iterator != 0){
		return iterator->get(idget, r2get, rget);
	}
}

int NeighborList::getLength(){
	return length;
}

void NeighborList::remove(int id){
	if(first == 0){
		return;
	}
	if(first->get() == id){
		Node* tmp = first;
		first = first->next;
		delete tmp;
	}else{
		first->remove(id);
	}
	length += -1;
}

bool NeighborList::exists(int id){
	for(Node* node = first; node != 0; node = node->next){
		if(node->get() == id){
			return true;
		}
	}
	return false;
}

void NeighborList::save(){
	saved = iterator;
}

void NeighborList::resume(){
	iterator = saved;
}

//New
void NeighborList::reset(){
	Node* nn;
	for(Node* n = first; n != 0; n = nn){
		nn = n->next;
		delete n;
	}
}

////////////////
//NODE METHODS//
////////////////

NeighborList::Node::Node(int id, double r2, double *r){
	this->id = id;
	this->r2 = r2;
	this->r = r;
	next = 0;
}

//Doesnt work? can not delete itself?
void NeighborList::Node::deleteAll(){
	if(next != 0){
		next->deleteAll();
	}
	delete this;
}

void NeighborList::Node::add(int val){
	if(next == 0){
		next = new Node(val);
	}else{
		next->add(val);
	}
}

void NeighborList::Node::remove(int remval){
	if(next == 0){
		return;
	}
	if(next->get() == remval){
		Node* tmp = next;
		next = next->next;
		delete tmp;
	}else{
		next->remove(remval);
	}
}

int NeighborList::Node::get(){
	return id;
}

void NeighborList::Node::get(int *idget, double *r2get, double *rget){
	*idget = id;
	*r2get = r2;
	rget = r;
}


