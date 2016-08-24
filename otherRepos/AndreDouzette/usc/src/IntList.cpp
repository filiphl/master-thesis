#include "IntList.h"

//TO DO: Implement proper iterator

IntList::IntList(){
	iterator = 0;
	first = 0;
	length = 0;
}

IntList::~IntList(){
	Node* nn;
	for(Node* n = first; n != 0; n = nn){
		nn = n->next;
		delete n;
	}
}

void IntList::add(int newval){
	if(first == 0){
		first = new Node(newval);
	}else{
		Node* tmp = first;
		first = new Node(newval);
		first->next = tmp;
	}
	length += 1;
}

int IntList::getFirst(){
	iterator = first;
	if(first != 0){
		return first->get();
	}else{
		return -1;
	}
}

bool IntList::hasNext(){
	return iterator != 0;
}

int IntList::getNext(){
	iterator = iterator->next;
	if(iterator != 0){
		return iterator->get();
	}else{
		return -1;
	}
}

int IntList::getLength(){
	return length;
}

void IntList::remove(int remval){
	if(first == 0){
		return;
	}
	if(first->get() == remval){
		Node* tmp = first;
		first = first->next;
		delete tmp;
	}else{
		first->remove(remval);
	}
	length += -1;
}

bool IntList::exists(int test){
	for(Node* node = first; node != 0; node = node->next){
		if(node->get() == test){
			return true;
		}
	}
	return false;
}

void IntList::save(){
	saved = iterator;
}

void IntList::resume(){
	iterator = saved;
}

//New
void IntList::reset(){
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

IntList::Node::Node(int val){
	this->val = val;
	next = 0;
}

//Doesnt work? can not delete itself?
void IntList::Node::deleteAll(){
	if(next != 0){
		next->deleteAll();
	}
	delete this;
}

void IntList::Node::add(int val){
	if(next == 0){
		next = new Node(val);
	}else{
		next->add(val);
	}
}

void IntList::Node::remove(int remval){
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

int IntList::Node::get(){
	return val;
}




