#include "MDcalculation.h"
using namespace std;

//Python wrapper for the c++ potential test program
int main(int argc, char** argv){
	//Test wether there are enough command line arguments
	if(argc < 2){
		//Insert not enough arguments exception?
		cout << argc << "\n";
		return 1;
	}
	
	//get path
	string path = argv[1];
	cout << path << "\n";
	string dataFilename = path + "/init/data.dat";
	MDcalculation* calculator = new MDcalculation(dataFilename, path, 0);
	calculator->potentialTest(path);
	calculator->finalize();
	return 0;
}