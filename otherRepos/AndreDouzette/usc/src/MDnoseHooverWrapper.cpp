#include "MDcalculation.h"
using namespace std;

//Python wrapper for the c++ simulation program
int main(int argc, char** argv){
	//Test wether there are enough command line arguments
	if(argc < 8){
		//Insert not enough arguments exception?
		cout << argc << "\n";
		return 1;
	}
		
	//////////////////////
	// Input parameters //
	//////////////////////
	//get path
	string path = argv[1];
	//get end index
	int length = atoi(argv[2]);
	//get save step
	int step = atoi(argv[3]);
	//get dt
	double dt = atof(argv[4]);
	//get tau
	double tau = atof(argv[5]);
	//get temp
	double T = atof(argv[6]);
	//get chain length
	int chainLength = atoi(argv[7]);
	
	length = length*step;
	/////////////////
	// Calculation //
	/////////////////
	MDcalculation* calculator = new MDcalculation(path, 0);
	// cout << "init complete\n";
	calculator->noseHoover(length, dt, tau, T, path, step, chainLength);
	calculator->finalize();
	
	return 0;
}