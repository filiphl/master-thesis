#ifndef INTLIST_H
#define INTLIST_H

class IntList{
	public:
		IntList();
		~IntList();
		void add(int newval);
		int getFirst();
		int getNext();
		bool hasNext();
		int getLength();
		void remove(int val);
		bool exists(int test);
		void save();
		void resume();
		void reset();
		
	private: class Node{
			private:
				int val;
			public:
				Node* next;
			public:
				Node(int val);
				void deleteAll();
				void add(int val);
				void remove(int val);
				int get();
		};
	private:
		Node* iterator;
		Node* first;
		int length;
		Node* saved;
};

#endif