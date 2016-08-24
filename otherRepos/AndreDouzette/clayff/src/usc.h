#ifndef USC_H
#define USC_H

//Need to import math for sqrt and exp

class usc{
	private:
		//2-body parameters
		bool* reacting2;
		double* H;
		double* W;
		double* D;
		double* Z;
		double* r1s;
		double* r4s;
		double rc;
		
		//3-body parameters
		bool* reacting3;
		double* r0;
		double* c0;
		double* xi;
		// double* C; This seems to be 0 in all situations
		double* B;
	public:
		
};

#endif