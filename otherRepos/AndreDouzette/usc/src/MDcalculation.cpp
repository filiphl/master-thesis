#include "MDcalculation.h"

using namespace std;

//!TO DO: Add option for silent mode?

//====================================================================//
// Misc.                                                              //
//====================================================================//

MDcalculation::MDcalculation(string path, int loadstate){
	Nchain = 0;
	//Initiate MPI, with some trash variables
	MPI_Init(new int[1], new char**[1]);
	//Get MPI id
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//Get MPI size
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	//Potential including scaling factors and force evaluation
	// cout << rank << "/" << procs << " Initialized\n";
	potential = new Potential();
	L0 = potential->L0;
	v0 = potential->v0;
	t0 = potential->t0;
	T0 = potential->T0;
	rc2 = potential->rc2;
	typeO = potential->typeO;
	typeH = potential->typeH;
	typeSi = potential->typeSi;
	typeOs = potential->typeOs;
	typeOh = potential->typeOh;
	dot0.reserve(100);
	dot1.reserve(100);
	ind0.reserve(100);
	ind1.reserve(100);
	//Read datafile, getting number of processors in each dimension and simulation space
	MPI_Barrier(MPI_COMM_WORLD);
	readDatafile(path + "/init/data.dat");
	MPI_Barrier(MPI_COMM_WORLD);
	readFile(path, loadstate);
	MPI_Barrier(MPI_COMM_WORLD);
	readMarkedFile(path + "/init/frozen.dat");
	readMonitoredFile(path + "/init/monitored.dat");
	// cout << rank << "/" << procs << " Starting simulation\n";
}

void MDcalculation::createDataStructures(){
	container = new ParticleContainer(Nparticles);
	//Number of particles by type
	NparticlesType = new int[3];
	//Mean displacement
	rMeanSquared = new double[3];
	rMeanSquaredSum = new double[3];
	
	double lx = sx/procsx;
	double ly = sy/procsy;
	double lz = sz/procsz;
	
	Nx = (int)(lx/potential->rc);
	Ny = (int)(ly/potential->rc);
	Nz = (int)(lz/potential->rc);
	
	Lx = lx/Nx;
	Ly = ly/Ny;
	Lz = lz/Nz;
	
	//Number of total unit cells
	// int cx = (int)(sx/potential->rc);
	// int cy = (int)(sy/potential->rc);
	// int cz = (int)(sz/potential->rc);
	// //Size of unit cells. L > rcut
	// Lx = sx/cx;
	// Ly = sy/cy;
	// Lz = sz/cz;
	// //Number of local cells
	// Nx = cx/procsx;
	// Ny = cy/procsy;
	// Nz = cz/procsz;
	//Processor coordinates
	zrank = rank/(procsx*procsy);
	yrank = rank/procsx - procsy*zrank;
	xrank = rank - procsx*(yrank + procsy*zrank);
	//Find neighbouring processors
	neighbour = new int[6];
	neighbour[0] = coordsToRank(xrank + 1, yrank, zrank);
	neighbour[1] = coordsToRank(xrank - 1, yrank, zrank);
	neighbour[2] = coordsToRank(xrank, yrank + 1, zrank);
	neighbour[3] = coordsToRank(xrank, yrank - 1, zrank);
	neighbour[4] = coordsToRank(xrank, yrank, zrank + 1);
	neighbour[5] = coordsToRank(xrank, yrank, zrank - 1);
	
	// //Adding extra cells to some processors
	// if(xrank < cx%procsx){
	// 	Nx++;
	// }
	// if(yrank < cy%procsy){
	// 	Ny++;
	// }
	// if(zrank < cz%procsz){
	// 	Nz++;
	// }
	
	//Create list
	//!To do: continous in memory, is it necessary?
	cells = new ParticleList***[Nx + 2];
	for(int i = 0; i < Nx + 2; i++){
		cells[i] = new ParticleList**[Ny + 2];
		for(int j = 0; j < Ny + 2; j++){
			cells[i][j] = new ParticleList*[Nz + 2];
			for(int k = 0; k < Nz + 2; k++){
				cells[i][j][k] = new ParticleList(container);
			}
		}
	}
	//Create particle export lists
	exportList = new ParticleList*[6];
	for(int i = 0; i < 6; i++){
		exportList[i] = new ParticleList(container);
	}
	//Create comunication lists
	communicateRight = new ParticleList(container);
	communicateLeft = new ParticleList(container);
	moveList = new ParticleList(container);
	b3bf2 = potential->r02[0][1]; //Largest 3BF cutoff
	b3bf = new double[4];
	b3bf[0] = -potential->r0[0][1];
	b3bf[1] = Lx*Nx - b3bf[0];
	b3bf[2] = Ly*Ny - b3bf[0];
	b3bf[3] = Lz*Nz - b3bf[0];
}

void MDcalculation::deleteLists(){
	for(int i = 0; i < Nx + 2; i++){
		for(int j = 0; j < Ny + 2; j++){
			for(int k = 0; k < Nz + 2; k++){
				// cells[i][j][k]->deleteParticles();
				delete cells[i][j][k];
			}
			delete[] cells[i][j];
		}
		delete[] cells[i];
	}
	delete[] cells;
	delete[] neighbour;
	for(int i = 0; i < 6; i++){
		// exportList[i]->deleteParticles();
		delete exportList[i];
	}
	delete[] exportList;
	delete communicateRight;
	delete communicateLeft;
	container->delAll();
	delete container;
	delete potential;
}

void MDcalculation::convertTime(double seconds, char* out){
	int s = (int)seconds;
	int m = (int)(s/60);
	int h = (int)(m/60);
	s = s - 60*m;
	m = m - 60*h;
	sprintf(out, "%02dh %02dm %02ds", h, m, s);
}

//====================================================================//
// File handling                                                      //
//====================================================================//

void MDcalculation::readDatafile(string filename){
	ifstream datafile;
	datafile.open(filename.c_str());
	if(datafile.is_open()){
		int t1, t2, t3;
		string line;
		char* strlist;
		//read processors px py pz
		getline(datafile, line);
		//split line
		strlist = strtok(strdup(line.c_str()), " ");
		procsx = stoi(strlist);
		strlist = strtok(NULL, " ");
		procsy = stoi(strlist);
		strlist = strtok(NULL, " ");
		procsz = stoi(strlist);
		//read sx sy sz
		getline(datafile, line);
		//split line
		strlist = strtok(strdup(line.c_str()), " ");
		sx = stod(strlist);
		strlist = strtok(NULL, " ");
		sy = stod(strlist);
		strlist = strtok(NULL, " ");
		sz = stod(strlist);
		datafile.close();
	}else{
		cout << "ERROR! CAN NOT OPEN DATA FILE!\n";
		//ADD ABORT
	}
}

void MDcalculation::readMarkedFile(string filename){
	ifstream file;
	string* linesplit;
	int nline;
	int id;
	file.open(filename.c_str());
	NnotFrozen = Nparticles;
	//Datafile does not exists means every particle is unmarked
	if(file.is_open()){
		string line;
		char* strlist;
		while(getline(file, line)){
			stringToArray(line, &linesplit, &nline);
			id = stoi(linesplit[0]);
			container->mark[id] = true;
			if(container->particles[id] != 0){
				container->particles[id]->v[0] = stod(linesplit[1])/v0;
				container->particles[id]->v[1] = stod(linesplit[2])/v0;
				container->particles[id]->v[2] = stod(linesplit[3])/v0;
			}
			NnotFrozen += -1;
		}
		file.close();
	}
}

void MDcalculation::readMonitoredFile(string filename){
	monitoring = false;
	ifstream file;
	string* linesplit;
	int nline;
	int id;
	file.open(filename.c_str());
	//Datafile does not exists means no particles are monitored
	if(file.is_open()){
		string line;
		char* strlist;
		while(getline(file, line)){
			stringToArray(line, &linesplit, &nline);
			id = stoi(linesplit[0]);
			container->monitored[id] = true;
			container->monitoredParticles.push_back(id);
			monitoring = true;
		}
		file.close();
	}
	if(monitoring){
		//Number of bins for each processor
		prN = 25;
		//Number of bins across all processors
		globalprN = prN*procsy;
		pressure = new double[globalprN];
		for(int i = 0; i < globalprN; i++){
			pressure[i] = 0;
		}
		//Length of bin
		prL = sy/(procsy*prN);
		//Inverse length
		prLi = 1/prL;
		pressure0 = potential->E0/(3*prL*sx*sz*potential->L0*potential->L0*potential->L0);
		//Start bin of this processor
		cy0 = prN*yrank;
		//Number of samples
		pressureSamples = 0;
	}
}

void MDcalculation::saveEnergy(string runName, int n){
	int pres = 10;
	ofstream energyFile;
	string name = runName + "/";
	char* nstr = new char[16];
	sprintf(nstr, "%05d/energy.dat", n);
	name.append(nstr);
	delete[] nstr;
	energyFile.open(name.c_str());
	if(energyFile.is_open()){
		energyFile << setprecision(pres) << K*potential->E0 << "\n";
		energyFile << setprecision(pres) << V*potential->E0 << "\n";
		energyFile << setprecision(pres) << 2*T0*K/(3.0*NnotFrozen) << "\n";
		energyFile << setprecision(pres) << rMeanSquaredSum[0]*L0*L0 << " ";
		energyFile << setprecision(pres) << rMeanSquaredSum[1]*L0*L0 << " ";
		energyFile << setprecision(pres) << rMeanSquaredSum[2]*L0*L0 << "\n";
		energyFile.close();
		for(int i = 0; i < Nchain; i++){
			energyFile << setprecision(pres) << zeta[i]/potential->t0 << " ";
		}
	}
}

void MDcalculation::saveInteratomicForces(string runName, int n){
	int pres = 15;
	ofstream forceFile;
	string name = runName + "/";
	char* nstr = new char[31];
	sprintf(nstr, "%05d/pressure_%03d.dat", n, rank);
	name.append(nstr);
	forceFile.open(name.c_str());
	delete[] nstr;
	int n1, n2;
	if(forceFile.is_open()){
		for(int i = 0; i < prN*procsy; i++){
			forceFile << pressure[i]*pressure0/pressureSamples << "\n";
		}
		forceFile.close();
	}
}

//TO DO: add error message if file format is wrong, and not able to read file
void MDcalculation::readFile(string filename, int n){
	string name = filename + "/";
	char* nstr = new char[13];
	sprintf(nstr, "%05d/%03d.xyz", n, rank);
	name.append(nstr);
	delete[] nstr;
	int tst = 0;
	
	double x, y, z;
	double vx, vy, vz;
	int cx, cy, cz;
	char* strlist;
	int i = 0;
	int id, type;
	string line;
	string* linesplit;
	int nline;
	
	int* Nlocal = new int[3];
	Nlocal[0] = 0; Nlocal[1] = 0; Nlocal[2] = 0;
	
	Particle* particle;
	
	ifstream stateFile;
	stateFile.open(name.c_str());
	MPI_Barrier(MPI_COMM_WORLD);
	if(stateFile.is_open()){
		//Read first line (number of particles)
		getline(stateFile, line);
		int N = atoi(line.c_str());
		//Gather number of particles from all processors
		MPI_Allreduce(&N, &Nparticles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		//Creating all datastructures for particle storrage, and processor neighbour lists
		createDataStructures();
		//read second line (comment)
		getline(stateFile, line);
		while(getline(stateFile, line)){
			stringToArray(line, &linesplit, &nline);
			//splits line by spaces
			//get type
			type = stoi(linesplit[0]);
			Nlocal[type] += 1;
			// //get position
			x = stod(linesplit[1])/L0;
			y = stod(linesplit[2])/L0;
			z = stod(linesplit[3])/L0;
			// if(x == 0) x += 0.0000001;
			// if(y == 0) y += 0.0000001;
			// if(z == 0) z += 0.0000001;
			// //Get velocity
			vx = stod(linesplit[4])/v0;
			vy = stod(linesplit[5])/v0;
			vz = stod(linesplit[6])/v0;
			//get id
			id = stoi(linesplit[7]);
			//create particle
			container->types[id] = type;
			particle = container->createParticle(x, y, z, vx, vy, vz, id);
			particle->rReal[0] = 0;
			particle->rReal[1] = 0;
			particle->rReal[2] = 0;
			//Add to cell, or export if out of bounds
			// if(x < 0){
			// 	exportList[1]->add(particle);
			// 	cout << "x: " << x << "\n";
			// }else if(x > Nx*Lx){
			// 	exportList[0]->add(particle);
			// 	cout << "X: " << x - Nx*Lx << "\n";
			// }else if(y < 0){
			// 	exportList[3]->add(particle);
			// 	cout << "y: " << y << "\n";
			// }else if(y > Ny*Ly){
			// 	exportList[2]->add(particle);
			// 	cout << "Y: " << y - Ny*Ly << "\n";
			// }else if(z < 0){
			// 	cout << "z: " << z << "\n";
			// 	exportList[5]->add(particle);
			// }else if(z > Nz*Lz){
			// 	exportList[4]->add(particle);
			// 	cout << "Z: " << z - Nz*Lz << "\n";
			// }else{
			
			//Adding to correct list.
			cx = (int)(x/Lx + 1);
			cy = (int)(y/Ly + 1);
			cz = (int)(z/Lz + 1);
			//Nudging particles on the boundary inside domain
			if(cx == Nx + 1){
				particle->r[0] += -0.00001;
				cx = Nx;
			}
			if(cy == Ny + 1){
				particle->r[1] += -0.00001;
				cy = Ny;
			}
			if(cz == Nz + 1){
				particle->r[2] += -0.00001;
				cz = Nz;
			}
			if(cx == 0){
				particle->r[0] += 0.00001;
				cx = 0;
			}
			if(cy == 0){
				particle->r[1] += 0.00001;
				cy = 0;
			}
			if(cz == 0){
				particle->r[2] += 0.00001;
				cz = 0;
			}
			
			cells[cx][cy][cz]->add(particle);
			
			// }
			//Adding particle to container
			delete[] linesplit;
		}
		stateFile.close();
		interchangeParticles();
	}else{
		cout << "ERROR! CAN NOT OPEN FILE!\n";
		//!ADD ABORT HERE
	}
	MPI_Allreduce(Nlocal, NparticlesType, 3, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	//Broadcasting the types of particles
	int* tmp = new int[Nparticles];
	MPI_Allreduce(container->types, tmp, Nparticles, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	delete[] container->types;
	container->types = tmp;
	delete[] Nlocal;
}

void MDcalculation::saveFile(string filename, int n){
	//To do: check if directory exists and create it if it doesnt
	string comment = "<comment>";
	int Np = 0;
	int pres = 10;
	ofstream stateFile;
	string name = filename + "/";
	char* nstr = new char[13];
	sprintf(nstr, "%05d/%03d.xyz", n, rank);
	name.append(nstr);
	delete[] nstr;
	stateFile.open(name.c_str());
	
	//Forces in eV/Aa
	double F0 = potential->E0/potential->L0;
	
	if(stateFile.is_open()){
		//get number of particles
		for(int i = 1; i < Nx + 1; i++){
		for(int j = 1; j < Ny + 1; j++){
		for(int k = 1; k < Nz + 1; k++){
			Np += cells[i][j][k]->getLength();
		}
		}
		}
		stateFile << Np << "\n" << comment.c_str() << "\n";
		for(int i = 1; i < Nx + 1; i++){
		for(int j = 1; j < Ny + 1; j++){
		for(int k = 1; k < Nz + 1; k++){
			for(Particle* p = cells[i][j][k]->getFirst(); cells[i][j][k]->hasNext(); p = cells[i][j][k]->getNext()){
				//particle type, change to chemical symbol?
				stateFile << p->type << " ";
				//position
				stateFile << setprecision(pres) << p->r[0]*L0;
				stateFile << " ";
				stateFile << setprecision(pres) << p->r[1]*L0;
				stateFile << " ";
				stateFile << setprecision(pres) << p->r[2]*L0;
				stateFile << " ";
				//velocity
				stateFile << setprecision(pres) << p->v[0]*v0;
				stateFile << " ";
				stateFile << setprecision(pres) << p->v[1]*v0;
				stateFile << " ";
				stateFile << setprecision(pres) << p->v[2]*v0;
				stateFile << " ";
				//particle id
				stateFile << p->id;
				stateFile << " ";
				stateFile << setprecision(pres) << p->F[0]*F0;
				stateFile << " ";
				stateFile << setprecision(pres) << p->F[1]*F0;
				stateFile << " ";
				stateFile << setprecision(pres) << p->F[2]*F0;
				//new line for each particle
				stateFile << "\n";
			}
		}
		}
		}
		stateFile.close();
	}else{
		cout << "ERROR: Can not open " << name << "\n";
	}
}

void MDcalculation::stringToArray(string str, string** strarray, int* length){
	stringstream ss;
	char* tok = strtok(strdup(str.c_str()), " ");
	*length = 0;
	while(tok != NULL){
		*length += 1;
		tok = strtok(NULL, " ");
	}
	*strarray = new string[*length];
	tok = strtok(strdup(str.c_str()), " ");
	for(int i = 0; i < *length; i++){
		(*strarray)[i] = *(new string(tok));
		tok = strtok(NULL, " ");
	}
	delete[] tok;
}

//====================================================================//
// Paralellization                                                    //
//====================================================================//

void MDcalculation::finalize(){
	//Temporary testing stuff
	
	//End of testing
	
	//To do: Delete lots of stuff here!
	
	//End of parallellization
	deleteLists();
	MPI_Finalize();
}

void MDcalculation::interchangeParticles(){
	//Deleting particles on boundary
	//This needs to be done first, such that no newly
	//imported particles are accidently deleted
	for(int j = 1; j < Ny + 1; j++){
		for(int k = 1; k < Nz + 1; k++){
			cells[Nx + 1][j][k]->deleteParticles();
			cells[0][j][k]->deleteParticles();
		}
	}
	for(int i = 0; i < Nx + 2; i++){
		for(int k = 1; k < Nz + 1; k++){
			cells[i][Ny + 1][k]->deleteParticles();
			cells[i][0][k]->deleteParticles();
		}
	}
	for(int i = 0; i < Nx + 2; i++){
		for(int j = 0; j < Ny + 2; j++){
			cells[i][j][Nz + 1]->deleteParticles();
			cells[i][j][0]->deleteParticles();
		}
	}
	//Interchanging particles in each dimension
	for(int dim = 0; dim < 3; dim++){
		interchangeParticles(dim);
	}
}

void MDcalculation::interchangeParticles(int dim){
	MPI_Status status;
	MPI_Request req0, req1;
	int* stats;
	int* statsIn;
	double* data;
	double* dataIn;
	int nStats, nData;
	int nStatsIn, nDataIn;
	//Senders and recievers
	int right = neighbour[2*dim];
	int left = neighbour[2*dim + 1];
	//Export arrays
	ParticleList* rightExport = exportList[2*dim];
	ParticleList* leftExport = exportList[2*dim + 1];
	
	//Communicating to right
	rightExport->makeArrays(&stats, &data, &nStats, &nData);
	rightExport->deleteParticles();
	shiftCoordinatesOut(data, nData, dim);
	//Send data to right
	MPI_Isend(data, nData, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &req0);
	//Recieve data from left
	MPI_Probe(left, 0, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_DOUBLE, &nDataIn);
	dataIn = new double[nDataIn];
	MPI_Irecv(dataIn, nDataIn, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &req1);
	MPI_Wait(&req0, &status);
	MPI_Wait(&req1, &status);
	//Send stats to right
	MPI_Isend(stats, nStats, MPI_INT, right, 1, MPI_COMM_WORLD, &req0);
	//Recieve stats from left
	MPI_Probe(left, 1, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_INT, &nStatsIn);
	statsIn = new int[nStatsIn];
	MPI_Irecv(statsIn, nStatsIn, MPI_INT, left, 1, MPI_COMM_WORLD, &req1);
	MPI_Wait(&req0, &status);
	MPI_Wait(&req1, &status);
	addParticles(statsIn, dataIn, nStatsIn, nDataIn);
	delete[] stats; delete[] data; delete[] statsIn; delete[] dataIn;
	
	//Communicating to left
	leftExport->makeArrays(&stats, &data, &nStats, &nData);
	leftExport->deleteParticles();
	//Send data to left
	MPI_Isend(data, nData, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &req0);
	//Recieve data from right
	MPI_Probe(right, 0, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_DOUBLE, &nDataIn);
	dataIn = new double[nDataIn];
	MPI_Irecv(dataIn, nDataIn, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &req1);
	MPI_Wait(&req0, &status);
	MPI_Wait(&req1, &status);
	//Send stats to left
	MPI_Isend(stats, nStats, MPI_INT, left, 1, MPI_COMM_WORLD, &req0);
	//Recieve stats from right
	MPI_Probe(right, 1, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_INT, &nStatsIn);
	statsIn = new int[nStatsIn];
	MPI_Irecv(statsIn, nStatsIn, MPI_INT, right, 1, MPI_COMM_WORLD, &req1);
	MPI_Wait(&req0, &status);
	MPI_Wait(&req1, &status);
	shiftCoordinatesIn(dataIn, nDataIn, dim);
	addParticles(statsIn, dataIn, nStatsIn, nDataIn);
	delete[] stats; delete[] data; delete[] statsIn; delete[] dataIn;
}

void MDcalculation::shiftCoordinatesOut(double* data, int nData, int dim){
	double L;
	if(dim == 0){
		L = Lx*Nx;
	}else if(dim == 1){
		L = Ly*Ny;
	}else{
		L = Lz*Nz;
	}
	for(int i = 0; i < nData/9; i++){
		data[i*9 + dim] += -L;
	}
}

void MDcalculation::shiftCoordinatesIn(double* data, int nData, int dim){
	double L;
	if(dim == 0){
		L = Lx*Nx;
	}else if(dim == 1){
		L = Ly*Ny;
	}else{
		L = Lz*Nz;
	}
	for(int i = 0; i < nData/9; i++){
		data[i*9 + dim] += L;
	}
}

void MDcalculation::shiftBoundaryCoordinatesOut(double* data, int nData, int dim){
	double L;
	if(dim == 0){
		L = Lx*Nx;
	}else if(dim == 1){
		L = Ly*Ny;
	}else{
		L = Lz*Nz;
	}
	for(int i = 0; i < nData/3; i++){
		data[i*3 + dim] += -L;
	}
}

void MDcalculation::shiftBoundaryCoordinatesIn(double* data, int nData, int dim){
	double L;
	if(dim == 0){
		L = Lx*Nx;
	}else if(dim == 1){
		L = Ly*Ny;
	}else{
		L = Lz*Nz;
	}
	for(int i = 0; i < nData/3; i++){
		data[i*3 + dim] += L;
	}
}

void MDcalculation::communicateBorder(){
	//Gathering lists in x direction
	for(int j = 1; j < Ny + 1; j++){
		for(int k = 1; k < Nz + 1; k++){
			communicateRight->add(cells[Nx][j][k]);
			communicateLeft->add(cells[1][j][k]);
		}
	}
	//Sending particles in x direction
	exportBoundary(communicateRight, communicateLeft, 0);
	//Reset lists
	communicateRight->reset();
	communicateLeft->reset();
	
	//Gathering lists in y direction
	for(int i = 0; i < Nx + 2; i++){
		for(int k = 1; k < Nz + 1; k++){
			communicateRight->add(cells[i][Ny][k]);
			communicateLeft->add(cells[i][1][k]);
		}
	}
	//Sending particles in y direction
	exportBoundary(communicateRight, communicateLeft, 1);
	//Reset lists
	communicateRight->reset();
	communicateLeft->reset();
	
	//Gathering lists in z direction
	for(int i = 0; i < Nx + 2; i++){
		for(int j = 0; j < Ny + 2; j++){
			communicateRight->add(cells[i][j][Nz]);
			communicateLeft->add(cells[i][j][1]);
		}
	}
	//Sending particles in z direction
	exportBoundary(communicateRight, communicateLeft, 2);
	//Reset lists
	communicateRight->reset();
	communicateLeft->reset();
}

void MDcalculation::exportBoundary(ParticleList* rightExport, ParticleList* leftExport, int dim){
	MPI_Status status;
	MPI_Request req0, req1;
	int right = neighbour[2*dim];
	int left  = neighbour[2*dim + 1];
	int* stats;
	double* data;
	int* statsIn;
	double* dataIn;
	int nData, nStats;
	int nDataIn, nStatsIn;
	
	//Communicating to right
	rightExport->makeBoundaryArrays(&stats, &data, &nStats, &nData);
	shiftBoundaryCoordinatesOut(data, nData, dim);
	//Send data to right
	MPI_Isend(data, nData, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &req0);
	//Recieve data from left
	MPI_Probe(left, 0, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_DOUBLE, &nDataIn);
	dataIn = new double[nDataIn];
	MPI_Irecv(dataIn, nDataIn, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &req1);
	MPI_Wait(&req0, &status);
	MPI_Wait(&req1, &status);
	//Send stats to right
	MPI_Isend(stats, nStats, MPI_INT, right, 1, MPI_COMM_WORLD, &req0);
	//Recieve stats from left
	MPI_Probe(left, 1, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_INT, &nStatsIn);
	statsIn = new int[nStatsIn];
	MPI_Irecv(statsIn, nStatsIn, MPI_INT, left, 1, MPI_COMM_WORLD, &req1);
	MPI_Wait(&req0, &status);
	MPI_Wait(&req1, &status);
	addBoundary(statsIn, dataIn, nStatsIn, nDataIn);
	delete[] stats; delete[] data; delete[] statsIn; delete[] dataIn;
	
	//Communicating to left
	leftExport->makeBoundaryArrays(&stats, &data, &nStats, &nData);
	//Send data to left
	MPI_Isend(data, nData, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &req0);
	//Recieve data from right
	MPI_Probe(right, 0, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_DOUBLE, &nDataIn);
	dataIn = new double[nDataIn];
	MPI_Irecv(dataIn, nDataIn, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &req1);
	MPI_Wait(&req0, &status);
	MPI_Wait(&req1, &status);
	//Send stats to left
	MPI_Isend(stats, nStats, MPI_INT, left, 1, MPI_COMM_WORLD, &req0);
	//Recieve stats from right
	MPI_Probe(right, 1, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_INT, &nStatsIn);
	statsIn = new int[nStatsIn];
	MPI_Irecv(statsIn, nStatsIn, MPI_INT, right, 1, MPI_COMM_WORLD, &req1);
	MPI_Wait(&req0, &status);
	MPI_Wait(&req1, &status);
	shiftBoundaryCoordinatesIn(dataIn, nDataIn, dim);
	addBoundary(statsIn, dataIn, nStatsIn, nDataIn);
	delete[] stats; delete[] data; delete[] statsIn; delete[] dataIn;
}

void MDcalculation::addBoundary(int* stats, double* data, int nStats, int nData){
	//Number of particles
	int N = nData/3;
	double x, y, z;
	int id, type;
	//Loop variable for stats array
	Particle* p;
	for(int i = 0; i < N; i++){
		//Find index for first entry of the current particle
		x = data[3*i];
		y = data[3*i + 1];
		z = data[3*i + 2];
		
		id = stats[i];
		//Add particle to cell list
		p = container->createParticle(x, y, z, id);
		addToBorder(p);
	}
}

void MDcalculation::addToBorder(Particle* p){
	//Adding particle to correct cell
	cells[(int)(p->r[0]/Lx + 1)][(int)(p->r[1]/Ly + 1)][(int)(p->r[2]/Lz + 1)]->add(p);
}

int MDcalculation::coordsToRank(int x, int y, int z){
	x = (x + procsx)%procsx;
	y = (y + procsy)%procsy;
	z = (z + procsz)%procsz;
	return x + procsx*(y + z*procsy);
}

void MDcalculation::addParticles(int* stats, double* data, int nStats, int nData){
	//Number of particles to add
	int N = nData/9;
	//Coordinates and velocity of new particles
	double x, y, z;
	double vx, vy, vz;
	//Information about the particle
	int id;
	//Loop variable for stats array
	Particle* p;
	for(int i = 0; i < N; i++){
		x = data[9*i];
		y = data[9*i + 1];
		z = data[9*i + 2];
		vx = data[9*i + 3];
		vy = data[9*i + 4];
		vz = data[9*i + 5];
		id = stats[i];
		p = container->createParticle(x, y, z, vx, vy, vz, id);
		p->rReal[0] = data[9*i + 6];
		p->rReal[1] = data[9*i + 7];
		p->rReal[2] = data[9*i + 8];
		addToSystem(p);
	}
}

void MDcalculation::addToSystem(Particle* p){
	double x = p->r[0];
	double y = p->r[1];
	double z = p->r[2];
	//Tests if particle is outside of simulation area,
	//then adds it to the correct cell list.
	if(x <= 0){
		exportList[1]->add(p);
	}else if(x > Nx*Lx){
		exportList[0]->add(p);
	}else if(y <= 0){
		exportList[3]->add(p);
	}else if(y > Ny*Ly){
		exportList[2]->add(p);
	}else if(z <= 0){
		exportList[5]->add(p);
	}else if(z > Nz*Lz){
		exportList[4]->add(p);
	}else{
		//Adding to correct list.
		//Need to shift by one due to boundary lists
		cells[(int)(x/Lx + 1)][(int)(y/Ly + 1)][(int)(z/Lz + 1)]->add(p);
	}
}

void MDcalculation::addOxygenToList(ParticleList* list, int i, int j, int k){
	ParticleList* cellList = cells[i][j][k];
	for(Particle* p = cellList->getFirst(); cellList->hasNext(); p = cellList->getNext()){
		if(p->type != typeO){
			continue;
		}
		list->add(p);
	}
}

void MDcalculation::communicateWeights(){
	//Gathering lists in x direction
	for(int j = 1; j < Ny + 1; j++){
		for(int k = 1; k < Nz + 1; k++){
			addOxygenToList(communicateRight, Nx, j, k);
			addOxygenToList(communicateLeft, 1, j, k);
		}
	}
	//Sending particles in x direction
	exportWeights(communicateRight, communicateLeft, 0);
	//Reset lists
	communicateRight->reset();
	communicateLeft->reset();
	
	//Gathering lists in y direction
	for(int i = 0; i < Nx + 2; i++){
		for(int k = 1; k < Nz + 1; k++){
			addOxygenToList(communicateRight, i, Ny, k);
			addOxygenToList(communicateLeft, i, 1, k);
		}
	}
	//Sending particles in y direction
	exportWeights(communicateRight, communicateLeft, 1);
	//Reset lists
	communicateRight->reset();
	communicateLeft->reset();
	
	//Gathering lists in z direction
	for(int i = 0; i < Nx + 2; i++){
		for(int j = 0; j < Ny + 2; j++){
			addOxygenToList(communicateRight, i, j, Nz);
			addOxygenToList(communicateLeft, i, j, 1);
		}
	}
	//Sending particles in z direction
	exportWeights(communicateRight, communicateLeft, 2);
	//Reset lists
	communicateRight->reset();
	communicateLeft->reset();
}

void MDcalculation::exportWeights(ParticleList* rightExport, ParticleList* leftExport, int dim){
	MPI_Status status;
	MPI_Request req0, req1;
	int right = neighbour[2*dim];
	int left  = neighbour[2*dim + 1];
	int* stats;		//ids
	double* data;	//nSi and nH
	int* statsIn;
	double* dataIn;
	int nData, nStats;
	int nDataIn, nStatsIn;
	
	//Communicating to right
	rightExport->makeWeightArrays(&stats, &data, &nStats, &nData, typeO);
	//Send data to right
	MPI_Isend(data, nData, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &req0);
	//Recieve data from left
	MPI_Probe(left, 0, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_DOUBLE, &nDataIn);
	dataIn = new double[nDataIn];
	MPI_Irecv(dataIn, nDataIn, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &req1);
	MPI_Wait(&req0, &status);
	MPI_Wait(&req1, &status);
	//Send stats to right
	MPI_Isend(stats, nStats, MPI_INT, right, 1, MPI_COMM_WORLD, &req0);
	//Recieve stats from left
	MPI_Probe(left, 1, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_INT, &nStatsIn);
	statsIn = new int[nStatsIn];
	MPI_Irecv(statsIn, nStatsIn, MPI_INT, left, 1, MPI_COMM_WORLD, &req1);
	MPI_Wait(&req0, &status);
	MPI_Wait(&req1, &status);
	addWeights(statsIn, dataIn, nStatsIn, nDataIn);
	delete[] stats; delete[] data; delete[] statsIn; delete[] dataIn;
	
	//Communicating to left
	leftExport->makeWeightArrays(&stats, &data, &nStats, &nData, typeO);
	//Send data to left
	MPI_Isend(data, nData, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &req0);
	//Recieve data from right
	MPI_Probe(right, 0, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_DOUBLE, &nDataIn);
	dataIn = new double[nDataIn];
	MPI_Irecv(dataIn, nDataIn, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &req1);
	MPI_Wait(&req0, &status);
	MPI_Wait(&req1, &status);
	//Send stats to left
	MPI_Isend(stats, nStats, MPI_INT, left, 1, MPI_COMM_WORLD, &req0);
	//Recieve stats from right
	MPI_Probe(right, 1, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_INT, &nStatsIn);
	statsIn = new int[nStatsIn];
	MPI_Irecv(statsIn, nStatsIn, MPI_INT, right, 1, MPI_COMM_WORLD, &req1);
	MPI_Wait(&req0, &status);
	MPI_Wait(&req1, &status);
	addWeights(statsIn, dataIn, nStatsIn, nDataIn);
	delete[] stats; delete[] data; delete[] statsIn; delete[] dataIn;
}

void MDcalculation::addWeights(int* stats, double* data, int nStats, int nData){
	for(int i = 0; i < nStats; i++){
		container->get(stats[i])->nSi += data[2*i];
		container->get(stats[i])->nH += data[2*i + 1];
	}
}

//====================================================================//
// Force calculation                                                  //
//====================================================================//

void MDcalculation::calculateForces(bool monitoring){
	V = 0;
	ParticleList* list;
	//Calculating weights, finding 2BF and 3BF particle lists
	findNeighbours();
	communicateWeights();
	//Calculating 2BF and 3BF
	if(monitoring){
		//Temp fix for monitoring particles.
		//It suddenly turns off for no reason! FIX!
		for(int i = 0; i < Nparticles; i++){
			if(container->particles[i] != 0 && !container->particles[i]->outside){
				container->particles[i]->monitored = container->monitored[i];
			}
		}
		for(int i = 0; i < Nx + 2; i++){
		for(int j = 0; j < Ny + 2; j++){
		for(int k = 0; k < Nz + 2; k++){
			list = cells[i][j][k];
			if(i != 0 && i != Nx + 1 && j != 0 && j != Ny + 1 && k != 0 && k != Nz + 1){
			//Internal 3BF calculation
				for(Particle* p = list->getFirst(); list->hasNext(); p = list->getNext()){
					V += potential->force2Monitoring(p);
					V += potential->force3(p);
				}
			}else{
				//3BF for boundary with internal
				for(Particle* p = list->getFirst(); list->hasNext(); p = list->getNext()){
					V += potential->boundaryForce3(p);
				}
			}
		}
		}
		}
	}else{
		for(int i = 0; i < Nx + 2; i++){
		for(int j = 0; j < Ny + 2; j++){
		for(int k = 0; k < Nz + 2; k++){
			list = cells[i][j][k];
			if(i != 0 && i != Nx + 1 && j != 0 && j != Ny + 1 && k != 0 && k != Nz + 1){
			//Internal 3BF calculation
				for(Particle* p = list->getFirst(); list->hasNext(); p = list->getNext()){
					V += potential->force2(p);
					V += potential->force3(p);
				}
			}else{
				//3BF for boundary with internal
				for(Particle* p = list->getFirst(); list->hasNext(); p = list->getNext()){
					V += potential->boundaryForce3(p);
				}
			}
		}
		}
		}
	}
}

void MDcalculation::addOxygenNeighbour(Particle* p, Particle* q){
	//p is oxygen, q is H or Si
	//Loop unrolling to optimize calculation
	dr[0] = q->r[0] - p->r[0];
	dr[1] = q->r[1] - p->r[1];
	dr[2] = q->r[2] - p->r[2];
	double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
	if(r2 <= rc2){
		p->neighbours.push_back(q);
		if(q->type == typeSi){
			p->nSi += potential->numberInterpolate(p, q, r2);
		}else{
			p->nH += potential->numberInterpolate(p, q, r2);
		}
		if(r2 <= potential->r02[p->type][q->type]){
			p->bonds.push_back(q);
			q->bonds.push_back(p);
		}
	}
}

void MDcalculation::addNeighbour(Particle* p, Particle* q){
	//Loop unrolling to optimize calculation
	dr[0] = q->r[0] - p->r[0];
	dr[1] = q->r[1] - p->r[1];
	dr[2] = q->r[2] - p->r[2];
	double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
	if(r2 <= rc2){
		p->neighbours.push_back(q);
		if(r2 <= potential->r02[p->type][q->type]){
			p->bonds.push_back(q);
			q->bonds.push_back(p);
		}
	}
}

void MDcalculation::addCloseBond(Particle* p, Particle* q){
	//Loop unrolling to optimize calculation
	dr[0] = q->r[0] - p->r[0];
	dr[1] = q->r[1] - p->r[1];
	dr[2] = q->r[2] - p->r[2];
	if(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] <= potential->r02[p->type][q->type]){
		p->bonds.push_back(q);
		q->bonds.push_back(p);
	}
}

void MDcalculation::quickSort(Particle** array, int f, int l){
	int i, j;
	Particle* pivot;
	Particle* tmp;
	while(f < l){
		i = f - 1, j = l + 1;
		Particle* tmp;
		pivot = array[f];
		while(true){
			while(pivot->dot < array[--j]->dot);
			while(array[++i]->dot < pivot->dot);
			if(i < j){
				tmp = array[i];
				array[i] = array[j];
				array[j] = tmp;
			}else{
				break;
			}
		}
		quickSort(array, f, j);
		f = j + 1;
	}
}

void MDcalculation::sortParticles(Particle** data, int length){
	
	// quickSort(data, 0, length - 1);
	// return;
	
	// sort(&(data[0]), &(data[length - 1]), [&](Particle* a, Particle* b){return a->dot > b->dot;});
	
	int j;
	Particle* tmp;
	for(int i = 1; i < length; i++){
		j = i;
		while(j > 0 && data[j]->dot < data[j - 1]->dot){
			tmp = data[j];
			data[j] = data[j - 1];
			data[j - 1] = tmp;
			j--;
		}
	}
}

void MDcalculation::findNeighboursInternal(int i0, int j0, int k0, int i1, int j1, int k1){
	Particle** list0 = cells[i0][j0][k0]->array;
	Particle** list1 = cells[i1][j1][k1]->array;
	int length0 = cells[i0][j0][k0]->length;
	int length1 = cells[i1][j1][k1]->length;
	int di = i1 - i0;
	int dj = j1 - j0;
	int dk = k1 - k0;
	for(int i = 0; i < length0; i++){
		list0[i]->dot = di*list0[i]->r[0] + dj*list0[i]->r[1] + dk*list0[i]->r[2];
	}
	for(int i = 0; i < length1; i++){
		list1[i]->dot = di*list1[i]->r[0] + dj*list1[i]->r[1] + dk*list1[i]->r[2];
	}
	
	// sortParticles(list0, length0);
	sortParticles(list1, length1);
	
	double cutoff = (di*di + dj*dj + dk*dk)*rc2;
	double d;
	Particle* p;
	Particle* q;
	for(int i = 0; i < length0; i++){
		p = list0[i];
		if(p->type == typeO){
			for(int j = 0; j < length1; j++){
			q = list1[j];
				d = q->dot - p->dot;
				if(d*d < cutoff){
					if(q->type == typeO){
						addNeighbour(p, q);
					}else{
						addOxygenNeighbour(p, q);
					}
				}else{
					break;
				}
			}
		}else{
			for(int j = 0; j < length1; j++){
			q = list1[j];
				d = q->dot - p->dot;
				if(d*d < cutoff){
					if(q->type == typeO){
						addOxygenNeighbour(q, p);
					}else{
						addNeighbour(p, q);
					}
				}else{
					break;
				}
			}
		}
	}
}

void MDcalculation::findNeighboursBetween(int i0, int j0, int k0, int i1, int j1, int k1){
	Particle** list0 = cells[i0][j0][k0]->array;
	Particle** list1 = cells[i1][j1][k1]->array;
	int length0 = cells[i0][j0][k0]->length;
	int length1 = cells[i1][j1][k1]->length;
	int di = i1 - i0;
	int dj = j1 - j0;
	int dk = k1 - k0;
	for(int i = 0; i < length0; i++){
		list0[i]->dot = di*list0[i]->r[0] + dj*list0[i]->r[1] + dk*list0[i]->r[2];
	}
	for(int i = 0; i < length1; i++){
		list1[i]->dot = di*list1[i]->r[0] + dj*list1[i]->r[1] + dk*list1[i]->r[2];
	}
	
	// sortParticles(list0, length0);
	sortParticles(list1, length1);
	
	double cutoff = (di*di + dj*dj + dk*dk)*rc2;
	double d;
	Particle* p;
	Particle* q;
	for(int i = 0; i < length0; i++){
		p = list0[i];
		if(p->type == typeO){
			for(int j = 0; j < length1; j++){
			q = list1[j];
				d = q->dot - p->dot;
				if(d*d < cutoff){
					if(q->type == typeO){
						addNeighbour(p, q);
					}else{
						addOxygenNeighbour(p, q);
					}
				}else{
					break;
				}
			}
		}else{
			for(int j = 0; j < length1; j++){
			q = list1[j];
				d = q->dot - p->dot;
				if(d*d < cutoff){
					addNeighbour(p, q);
				}else{
					break;
				}
			}
		}
	}
}

void MDcalculation::findNeighboursBoundary(int i0, int j0, int k0, int i1, int j1, int k1){
	Particle** list0 = cells[i0][j0][k0]->array;
	Particle** list1 = cells[i1][j1][k1]->array;
	int length0 = cells[i0][j0][k0]->length;
	int length1 = cells[i1][j1][k1]->length;
	if(length0 == 0 || length1 == 0) return;
	int di = i1 - i0;
	int dj = j1 - j0;
	int dk = k1 - k0;
	for(int i = 0; i < length0; i++){
		list0[i]->dot = di*list0[i]->r[0] + dj*list0[i]->r[1] + dk*list0[i]->r[2];
	}
	for(int i = 0; i < length1; i++){
		list1[i]->dot = di*list1[i]->r[0] + dj*list1[i]->r[1] + dk*list1[i]->r[2];
	}
	
	// sortParticles(list0, length0);
	sortParticles(list1, length1);
	
	double cutoff = (di*di + dj*dj + dk*dk)*b3bf2;
	double d;
	Particle* p;
	Particle* q;
	for(int i = 0; i < length0; i++){
		p = list0[i];
		for(int j = 0; j < length1; j++){
			q = list1[j];
			d = q->dot - p->dot;
			if(d*d < cutoff){
				addCloseBond(p, q);
			}else{
				break;
			}
		}
	}
}

void MDcalculation::findNeighbours(){
	ParticleList* list1;
	ParticleList* list2;
	int i2, j2, k2;
	//Resetting number of H and Si neighbours of oxygens and monitoring lists
	for(int i = 1; i < Nx + 1; i++){
	for(int j = 1; j < Ny + 1; j++){
	for(int k = 1; k < Nz + 1; k++){
		list1 = cells[i][j][k];
		for(Particle* p = list1->getFirst(); list1->hasNext(); p = list1->getNext()){
			if(p->type == typeO){
				p->nH = 0;
				p->nSi = 0;
			}
			//Deleting bonds and neighbour lists
			p->resetLists();
		}
	}
	}
	}
	for(int i = 0; i < Nx + 2; i++){
	for(int j = 0; j < Ny + 2; j++){
	for(int k = 0; k < Nz + 2; k++){
		cells[i][j][k]->makeArray();
	}
	}
	}
	
	//Internal cells
	for(int i = 1; i < Nx + 1; i++){
	for(int j = 1; j < Ny + 1; j++){
	for(int k = 1; k < Nz + 1; k++){
		list1 = cells[i][j][k];
		//Calculating internal neighbours
		for(Particle* p = list1->getFirst(); list1->hasNext(); p = list1->getNext()){
			if(p->type == typeO){
				for(Particle* q = list1->getFirst(); p->id != q->id; q = list1->getNext()){
					if(q->type != typeO){
						addOxygenNeighbour(p, q);
					}else{
						addNeighbour(p, q);
					}
				}
			}else{
				for(Particle* q = list1->getFirst(); p->id != q->id; q = list1->getNext()){
					if(q->type == typeO){
						addOxygenNeighbour(q, p);
					}else{
						addNeighbour(p, q);
					}
				}
			}
		}
		//Half of surrounding cells (N3L)
		for(int n = 14; n < 27; n++){
			i2 = i + n%3 - 1;
			j2 = j + (n/3)%3 - 1;
			k2 = k + n/9 - 1;
			if(i2 != 0 && i2 != Nx + 1 && j2 != 0 && j2 != Ny + 1 && k2 != 0 && k2 != Nz + 1){
				findNeighboursInternal(i, j, k, i2, j2, k2);
			}else{
				findNeighboursBetween(i, j, k, i2, j2, k2);
			}
		}
	}
	}
	}
	
	//Border cells
	for(int i = 0; i < Nx + 2; i++){
	for(int j = 0; j < Ny + 2; j++){
	for(int k = 0; k < Nz + 2; k++){
		//Internal cell
		if(i != 0 && i != Nx + 1 && j != 0 && j != Ny + 1 && k != 0 && k != Nz + 1){
			continue;
		}
		list1 = cells[i][j][k];
		//Internal cell neighbours
		for(Particle* p = list1->getFirst(); list1->hasNext(); p = list1->getNext()){
			for(Particle* q = list1->getFirst(); q->id != p->id; q = list1->getNext()){
				addCloseBond(p, q);
			}
		}
		for(int n = 14; n < 27; n++){
			i2 = i + n%3 - 1;
			j2 = j + (n/3)%3 - 1;
			k2 = k + n/9 - 1;
			//out of bounds
			if(i2 == -1 || i2 == Nx + 2 || j2 == -1 || j2 == Ny + 2 || k2 == -1 || k2 == Nz + 2){
				continue;
			}
			if(i2 != 0 && i2 != Nx + 1 && j2 != 0 && j2 != Ny + 1 && k2 != 0 && k2 != Nz + 1){
				findNeighboursBetween(i2, j2, k2, i, j, k);
			}else{
				findNeighboursBoundary(i, j, k, i2, j2, k2);
			}
			
		}
	}
	}	
	}
}

//====================================================================//
// Integration                                                        //
//====================================================================//

void MDcalculation::run(int Nt, double dt, string runName, int saveInterval){
	potential->createForceFactor(dt);
	//Communicate out of local space particles
	interchangeParticles();
	communicateBorder();
	double time0, time1;
	int pc;
	time0 = MPI_Wtime();
	char* timeleft = new char[11];
	
	//Converting to numerical values
	dt = dt/t0;
	//Initialize forces
	double Vsum;
	calculateForces(false);
	//Saving energy for init state
	calculateKineticEnergy();
	calculateMeanDisplacement();
	MPI_Allreduce(&V, &Vsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	V = Vsum;
	if(rank == 0){
		saveEnergy(runName, 0);
	}
	if(rank == 0){
		cout << "[sim]   0.00 complete (??h ??m ??s)\n";
	}
	//Computation loop
	for(int i = 1; i <= Nt; i++){
		if(i%saveInterval == 0){
			integrate(dt, monitoring);
			saveFile(runName, i/saveInterval);
			//Saving energy to file
			MPI_Allreduce(&V, &Vsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			V = Vsum;
			calculateKineticEnergy();
			calculateMeanDisplacement();
			calculatePressure();
			if(rank == 0){
				saveEnergy(runName, i/saveInterval);
			}
			if(monitoring){
				saveInteratomicForces(runName, i/saveInterval);
			}
		}else{
			integrate(dt, false);
			calculatePressure();
		}
		if(rank == 0){
			//Print each percents completeness, needs to change! Not working
			//for values which are not a multiple of 100!
			if(i%((int)(0.01*Nt)) == 0){
				time1 = MPI_Wtime();
				pc = (int)((100*i)/Nt);
				if(pc != 0){
					convertTime((1.2 - 0.2*pc*0.01)*(time1 - time0)*(100 - pc)/pc, timeleft);
					printf("\033[1A\033[K[sim]   %.2f complete (%s)\n", pc*0.01, timeleft);
				}
			}
		}
	}
	if(rank == 0){
		convertTime(MPI_Wtime() - time0, timeleft);
		printf("\033[1A\033[K[sim]   Elapsed time: %s\n", timeleft);
	}
	delete[] timeleft;
}

void MDcalculation::thermalize(int Nt, double dt, double tau, int avglength, double T, string runName, int saveInterval){
	dt = dt/t0;
	potential->createForceFactor(dt);
	//Communicate out of local space particles
	interchangeParticles();
	communicateBorder();
	double time0, time1;
	int pc;
	time0 = MPI_Wtime();
	char* timeleft = new char[11];
	
	//Converting to numerical values
	T = T/T0;
	//Initialize forces
	double Vsum;
	calculateForces(false);
	// return;
	//Saving energy for init state
	calculateKineticEnergy();
	//average array initialization
	double* Tv = new double[avglength];
	double Tavg = K/(1.5*NnotFrozen);
	for(int i = 0; i < avglength; i++){
		Tv[i] = Tavg;
	}
	calculateMeanDisplacement();
	MPI_Allreduce(&V, &Vsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	V = Vsum;
	if(rank == 0){
		saveEnergy(runName, 0);
	}
	//Computation loop
	timeForce = 0;
	timeInt = 0;
	timeCom = 0;
	timeOther = 0;
	// double tmptime;
	// double timeTotal;
	if(rank == 0){
		cout << "[therm] 0.00 complete (??h ??m ??s)\n";
	}
	for(int i = 1; i <= Nt; i++){
		// tmptime = MPI_Wtime();
		if(i%saveInterval == 0){
			integrate(dt, monitoring);
			saveFile(runName, i/saveInterval);
			//Saving energy to file
			MPI_Allreduce(&V, &Vsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			V = Vsum;
			calculateKineticEnergy();
			calculateMeanDisplacement();
			calculatePressure();
			if(rank == 0){
				saveEnergy(runName, i/saveInterval);
			}
			if(monitoring){
				saveInteratomicForces(runName, i/saveInterval);
			}
		}else{
			integrate(dt, false);
		}
		//Finding mean temperature
		calculateKineticEnergy();
		Tv[i%avglength] = K/(1.5*NnotFrozen);
		Tavg = 0;
		for(int j = 0; j < avglength; j++){
			Tavg += Tv[j];
		}
		Tavg = Tavg/avglength;
		//Thermalize system
		thermostat(dt, T, tau, Tavg);
		if(rank == 0){
			//Print each percents completeness, needs to change! Not working
			//for values which are not a multiple of 100!
			if(i%((int)(0.01*Nt)) == 0){
				time1 = MPI_Wtime();
				pc = (int)((100*i)/Nt);
				if(pc != 0){
					convertTime((1.2 - 0.2*pc*0.01)*(time1 - time0)*(100 - pc)/pc, timeleft);
					printf("\033[1A\033[K[therm] %.2f complete (%s)\n", pc*0.01, timeleft);
					// timeTotal = timeForce + timeCom + timeInt + timeOther;
					// cout << "Force calculation time: " << timeForce << "   \t" << timeForce/timeTotal << "\n";
					// cout << "Communication time:     " << timeCom << "   \t" << timeCom/timeTotal << "\n";
					// cout << "Integration time:       " << timeInt << "   \t" << timeInt/timeTotal << "\n";
					// cout << "Other:                  " << timeOther << "   \t" << timeOther/timeTotal << "\n";
					// // cout << "Max neighbours:         " << neighbours << "\n";
				}
			}
		}
	}
	if(rank == 0){
		convertTime(MPI_Wtime() - time0, timeleft);
		printf("\033[1A\033[K[therm] Elapsed time: %s\n", timeleft);
	}
	delete[] timeleft;
}

void MDcalculation::noseHoover(int Nt, double dt, double tau, double T, string runName, int saveInterval, int Nchain){
	dt = dt/t0;
	potential->createForceFactor(dt);
	double time0, time1;
	int pc;
	this->Nchain = Nchain;
	//Communicate out of local space particles
	interchangeParticles();
	communicateBorder();
	time0 = MPI_Wtime();
	char* timeleft = new char[11];
	
	zeta = new double[Nchain];
	mA = new double[Nchain];
	mB = new double[Nchain];
	mC = new double[Nchain];
	mD = new double[Nchain];
	
	for(int i = 0; i < Nchain; i++){
		zeta[i] = 0;
		mA[i] = 0;
		mB[i] = 0;
		mC[i] = 0;
		mD[i] = 0;
	}
	
	//Converting to numerical values
	T = T/T0;
	tau = tau/t0;
	//Degrees of freedom
	Qj = T*tau*tau;
	Q0 = Qj*NnotFrozen;
	Qji = 1/Qj;
	Q0i = 1/Q0;
	
	//Initialize forces
	double Vsum;
	calculateForces(false);
	//Saving energy for init state
	calculateKineticEnergy();
	calculateMeanDisplacement();
	MPI_Allreduce(&V, &Vsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	V = Vsum;
	if(rank == 0){
		saveEnergy(runName, 0);
	}
	if(rank == 0){
		cout << "[noseH] 0.00 complete (??h ??m ??s)\n";
	}
	//Computation loop
	for(int i = 1; i <= Nt; i++){
		if(i%saveInterval == 0){
			integrateNoseHoover(dt, T, monitoring);
			saveFile(runName, i/saveInterval);
			//Saving energy to file
			MPI_Allreduce(&V, &Vsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			V = Vsum;
			calculateKineticEnergy();
			calculateMeanDisplacement();
			calculatePressure();
			if(rank == 0){
				saveEnergy(runName, i/saveInterval);
			}
			if(monitoring){
				saveInteratomicForces(runName, i/saveInterval);
			}
		}else{
			integrateNoseHoover(dt, T, false);
		}
		if(rank == 0){
			//Print each percents completeness, needs to change! Not working
			//for values which are not a multiple of 100!
			if(i%((int)(0.01*Nt)) == 0){
				time1 = MPI_Wtime();
				pc = (int)((100*i)/Nt);
				if(pc != 0){
					convertTime((1.2 - 0.2*pc*0.01)*(time1 - time0)*(100 - pc)/pc, timeleft);
					printf("\033[1A\033[K[noseH] %.2f complete (%s)\n", pc*0.01, timeleft);
				}
			}
		}
	}
	if(rank == 0){
		convertTime(MPI_Wtime() - time0, timeleft);
		printf("\033[1A\033[K[noseH] Elapsed time: %s\n", timeleft);
	}
	delete[] timeleft;
}

int MDcalculation::globalParticles(){
	//Finds number of particles across all processors (Used for bugtesting)
	int pl, pg;
	pl = localParticles();
	MPI_Allreduce(&pl, &pg, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	return pg;
}

int MDcalculation::localParticles(){
	//Finds number of local processor particles (Used for bugtesting)
	int pl = 0;
	for(int i = 0; i < Nx + 1; i++){
	for(int j = 0; j < Ny + 1; j++){
	for(int k = 0; k < Nz + 1; k++){
		pl += cells[i][j][k]->getLength();
	}
	}
	}
	for(int i = 0; i < 6; i++){
		pl += exportList[i]->getLength();
	}
	return pl;
}

void MDcalculation::integrate(double dt, bool monitoring){
	int p1, p2;
	int cellx, celly, cellz;
	//reset export lists, is this necessary here?
	for(int d = 0; d < 6; d++){
		exportList[d]->deleteParticles();
	}
	//Calculate new positions
	for(int i = 1; i <= Nx; i++){
	for(int j = 1; j <= Ny; j++){
	for(int k = 1; k <= Nz; k++){
		for(Particle* p = cells[i][j][k]->getFirst(); cells[i][j][k]->hasNext(); p = cells[i][j][k]->getNext()){
			//Calculate half time velcotites
			if(!container->mark[p->id]){
				potential->updateV(p);
			}
			//Resetting forces before new force calculation
			p->resetForces();
			//Calculate new positions
			potential->updateR(p, dt);
			//Move particles to export lists if outside of local simulation space
			if(p->r[0] > Nx*Lx){
				exportList[0]->add(p);
				cells[i][j][k]->remove(p->id);
			}else if(p->r[0] < 0){
				exportList[1]->add(p);
				cells[i][j][k]->remove(p->id);
			}else if(p->r[1] > Ny*Ly){
				exportList[2]->add(p);
				cells[i][j][k]->remove(p->id);
			}else if(p->r[1] < 0){
				exportList[3]->add(p);
				cells[i][j][k]->remove(p->id);
			}else if(p->r[2] > Nz*Lz){
				exportList[4]->add(p);
				cells[i][j][k]->remove(p->id);
			}else if(p->r[2] < 0){
				exportList[5]->add(p);
				cells[i][j][k]->remove(p->id);
			}else{
				//Finds new cell
				cellx = (int)(p->r[0]/Lx + 1);
				celly = (int)(p->r[1]/Ly + 1);
				cellz = (int)(p->r[2]/Lz + 1);
				//Moves particle to new cell
				if(i != cellx || j != celly || k != cellz){
					moveList->add(p);
					cells[i][j][k]->remove(p->id);
				}
			}
		}
	}
	}
	}
	//Adding moved particles to new cells
	for(Particle* p = moveList->getFirst(); moveList->hasNext(); p = moveList->getNext()){
		cellx = (int)(p->r[0]/Lx + 1);
		celly = (int)(p->r[1]/Ly + 1);
		cellz = (int)(p->r[2]/Lz + 1);
		cells[cellx][celly][cellz]->add(p);
	}
	moveList->reset();
	//Communicate out of local space particles
	interchangeParticles();
	communicateBorder();
	calculateForces(monitoring);
	//Calculate final velocities
	for(int i = 1; i <= Nx; i++){
	for(int j = 1; j <= Ny; j++){
	for(int k = 1; k <= Nz; k++){
		for(Particle* p = cells[i][j][k]->getFirst(); cells[i][j][k]->hasNext(); p = cells[i][j][k]->getNext()){
			if(!container->mark[p->id]){
				potential->updateV(p);
			}
		}
	}
	}
	}
}

void MDcalculation::integrateNoseHoover(double dt, double T, bool monitoring){
	int p1, p2;
	int cellx, celly, cellz;
	double tmp;
	
	calculateKineticEnergy();
	//reset export lists, is this necessary here?
	for(int d = 0; d < 6; d++){
		exportList[d]->deleteParticles();
	}
	//Calculate new positions
	for(int i = 1; i <= Nx; i++){
	for(int j = 1; j <= Ny; j++){
	for(int k = 1; k <= Nz; k++){
		for(Particle* p = cells[i][j][k]->getFirst(); cells[i][j][k]->hasNext(); p = cells[i][j][k]->getNext()){
			//Calculate half time velcotites
			if(!container->mark[p->id]){
				potential->updateVnh1(p, zeta[0]);
			}
			//Resetting forces before new force calculation
			p->resetForces();
			//Calculate new positions
			potential->updateR(p, dt);
			
			//Move particles to export lists if outside of local simulation space
			if(p->r[0] > Nx*Lx){
				exportList[0]->add(p);
				cells[i][j][k]->remove(p->id);
			}else if(p->r[0] < 0){
				exportList[1]->add(p);
				cells[i][j][k]->remove(p->id);
			}else if(p->r[1] > Ny*Ly){
				exportList[2]->add(p);
				cells[i][j][k]->remove(p->id);
			}else if(p->r[1] < 0){
				exportList[3]->add(p);
				cells[i][j][k]->remove(p->id);
			}else if(p->r[2] > Nz*Lz){
				exportList[4]->add(p);
				cells[i][j][k]->remove(p->id);
			}else if(p->r[2] < 0){
				exportList[5]->add(p);
				cells[i][j][k]->remove(p->id);
			}else{
				//Finds new cell
				cellx = (int)(p->r[0]/Lx + 1);
				celly = (int)(p->r[1]/Ly + 1);
				cellz = (int)(p->r[2]/Lz + 1);
				//Moves particle to new cell
				if(i != cellx || j != celly || k != cellz){
					moveList->add(p);
					cells[i][j][k]->remove(p->id);
				}
			}
		}
	}
	}
	}
	//Adding moved particles to new cells
	for(Particle* p = moveList->getFirst(); moveList->hasNext(); p = moveList->getNext()){
		cellx = (int)(p->r[0]/Lx + 1);
		celly = (int)(p->r[1]/Ly + 1);
		cellz = (int)(p->r[2]/Lz + 1);
		cells[cellx][celly][cellz]->add(p);
	}
	moveList->reset();
	//Communicate out of local space particles
	interchangeParticles();
	communicateBorder();
	calculateKineticEnergy();
	//Calculating forces
	calculateForces(monitoring);
	//Nose-Hoover chain thermostat algorithm
	if(Nchain == 1){
		zeta[0] = zeta[0] + dt*Q0i*(K - 1.5*NnotFrozen*T);
	}else{
		mD[0] = zeta[0] - dt*Q0i*(1.5*NnotFrozen*T - K);
		mB[0] = 1 + 0.5*dt*zeta[1];
		mC[0] = 0.5*dt*zeta[0];
		for(int i = 1; i < Nchain; i++){
			if(i == 1){
				mA[i] = -0.5*dt*zeta[0]*Q0*Qji;
			}else{
				mA[i] = -0.5*dt*zeta[i - 1];
			}
			if(i != Nchain - 1){
				mB[i] = 1 + 0.5*dt*zeta[i + 1];
			}else{
				mB[i] = 1;
			}
			mC[i] = 0.5*dt*zeta[i];
			mD[i] = zeta[i] - 0.5*T*dt*Qji;
		}
		
		mC[0] = mC[0]/mB[0];
		mD[0] = mD[0]/mB[0];
		for(int i = 1; i < Nchain; i++){
			tmp = 1/(mB[i] - mA[i]*mC[i - 1]);
			mC[i] = mC[i]*tmp;
			mD[i] = (mD[i] - mA[i]*mD[i - 1])*tmp;
		}
		zeta[Nchain - 1] = mD[Nchain - 1];
		for(int i = Nchain - 2; i >= 0; i--){
			zeta[i] = mD[i] - mC[i]*zeta[i + 1];
		}
	}
	potential->createVelocityFactor(dt, zeta[0]);
	//Calculate final velocities
	for(int i = 1; i <= Nx; i++){
	for(int j = 1; j <= Ny; j++){
	for(int k = 1; k <= Nz; k++){
		for(Particle* p = cells[i][j][k]->getFirst(); cells[i][j][k]->hasNext(); p = cells[i][j][k]->getNext()){
			if(!container->mark[p->id]){
				potential->updateVnh2(p);
			}
		}
	}
	}
	}
}

void MDcalculation::thermostat(double dt, double T0, double tau, double T){
	double lambda = sqrt(1 + dt/tau*(T0/T - 1));
	//Scaling velocities
	ParticleList* list;
	for(int i = 1; i <= Nx; i++){
	for(int j = 1; j <= Ny; j++){
	for(int k = 1; k <= Nz; k++){
		list = cells[i][j][k];
		for(Particle* p = list->getFirst(); list->hasNext(); p = list->getNext()){
			if(container->mark[p->id]){
				continue;
			}
			for(int d = 0; d < 3; d++){
				p->v[d] = p->v[d]*lambda;
			}
		}
	}
	}
	}
}

//====================================================================//
// Analyzing                                                          //
//====================================================================//

void MDcalculation::calculateKineticEnergy(){
	K = 0;
	double Ksum;
	double v2;
	ParticleList* list;
	for(int i = 1; i < Nx + 1; i++){
	for(int j = 1; j < Ny + 1; j++){
	for(int k = 1; k < Nz + 1; k++){
		list = cells[i][j][k];
		for(Particle* p = list->getFirst(); list->hasNext(); p = list->getNext()){
			if(container->mark[p->id]){
				continue;
			}
			v2 = p->v[0]*p->v[0] + p->v[1]*p->v[1] + p->v[2]*p->v[2];
			K += potential->m[p->type]*v2;
		}
	}
	}
	}
	K = 0.5*K;
	MPI_Allreduce(&K, &Ksum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	K = Ksum;
}

void MDcalculation::calculateMeanDisplacement(){
	//Resetting
	for(int t = 0; t < 3; t++){
		rMeanSquared[t] = 0;
	}
	ParticleList* list;
	for(int i = 1; i < Nx + 1; i++){
	for(int j = 1; j < Ny + 1; j++){
	for(int k = 1; k < Nz + 1; k++){
		list = cells[i][j][k];
		for(Particle* p = list->getFirst(); list->hasNext(); p = list->getNext()){
			for(int d = 0; d < 3; d++){
				rMeanSquared[p->type] += p->rReal[d]*p->rReal[d];
			}
		}
	}
	}
	}
	MPI_Reduce(rMeanSquared, rMeanSquaredSum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(rank == 0){
		for(int t = 0; t < 3; t++){
			if(NparticlesType[t] == 0){
				rMeanSquaredSum[t] = 0;
			}else{
				rMeanSquaredSum[t] = rMeanSquaredSum[t]/NparticlesType[t];
			}
		}
	}
}

void MDcalculation::virial(double y1, double y2, double w){
	int Y1 = (int)(y1*prLi);
	int Y2 = (int)(y2*prLi);
	if(Y1 == Y2){
		pressure[Y1 + cy0] += w;
		return;
	}
	//y1 being the smallest
	if(Y1 > Y2){
		double tmp = y1;
		int tmpInt = Y1;
		y1 = y2;
		y2 = tmp;
		Y1 = Y2;
		Y2 = tmpInt;
	}
	double lineLengthInverseVirial = w/(y2 - y1);
	pressure[(Y1 + cy0 + globalprN)%globalprN] += (Y1*prL + 1 - y1)*lineLengthInverseVirial;
	pressure[(Y2 + cy0)%globalprN]             += (y2 - Y2*prL)*lineLengthInverseVirial;
	w = prL*lineLengthInverseVirial;
	for(int i = Y1 + 1; i < Y2; i++){
		pressure[(i + cy0 + globalprN)%globalprN] += w;
	}
}

void MDcalculation::calculatePressure(){
	if(!monitoring){
		return;
	}
	Particle* p;
	Particle* q;
	int qi;
	double w;
	pressureSamples += 1;
	for(int pi: container->monitoredParticles){
		p = container->particles[pi];
		if(p == 0 || p->outside){
			continue;
		}
		//Ideal gas contribution
		pressure[(int)(p->r[1]*prLi) + cy0] += potential->m[p->type]*(p->v[0]*p->v[0] + p->v[1]*p->v[1] + p->v[2]*p->v[2]);
		for(int i = 0; i < p->index2B.size(); i++){
			qi = p->index2B[i];
			//Only count particle pairs once
			if(container->monitored[qi] && pi > qi){
				continue;
			}
			q = container->particles[qi];
			//Virial contribution
			w = p->forces2B[3*i]*(q->r[0] - p->r[0]) + p->forces2B[3*i + 1]*(q->r[1] - q->r[1]) + p->forces2B[3*i + 2]*(q->r[2] - p->r[2]);
			virial(p->r[1], q->r[1], -w);
		}
	}
}



