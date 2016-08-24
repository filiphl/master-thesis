#include "MDcalculation.h"

using namespace std;

int main(int argc, char** argv){
	string dataFilename = "../cristobaliteTest/init/data.dat";
	MDcalculation* calculator = new MDcalculation(dataFilename);
	calculator->readFile("../cristobaliteTest", 0);
	calculator->run(800, 0.01, "../cristobaliteTest/", 1);
	calculator->finalize();
	return 0;
}