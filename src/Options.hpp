#ifndef _OPTIONSH_
#define _OPTIONSH_

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <boost/program_options.hpp>

#include "BC.hpp"

using namespace std;
namespace po = boost::program_options;

class Options{

    public:
	
	//Variable map...
	po::variables_map vm;

	//Root rank
	int root;

	//Time Stepping Information
	double CFL, maxTime;
	int maxTimeStep, filterStep, checkStep, dumpStep, imageStep;

	//Boundary Condition Stuff
	string bcXType_str, bcYType_str, bcZType_str;
	BC::BCType bcXType, bcYType, bcZType;

	string bcX0_str, bcX1_str, bcY0_str, bcY1_str, bcZ0_str, bcZ1_str;
	BC::BCKind bcX0, bcX1, bcY0, bcY1, bcZ0, bcZ1;
	double periodicDisp[3][3];

	//DOMAIN Information
	int Nx, Ny, Nz;

	//Pencil Decomposition Stuff
	int pRow, pCol;

	//SOLVER Information
	double alphaF, mu_ref;
	bool useTiming;

	//Restart stuff
	bool fromRestart;
	string filename;
	int startingTimestep;
    
    Options(int mpiRank){
   
      root = 0;
 
      if(mpiRank == root){
	ifstream input_file("MPI_3DCompact.in");

	po::options_description input("Input Parameters");

	input.add_options()
	    ("TIMESTEPPING.CFL", 	 po::value<double>(&CFL), "CFL Number")
	    ("TIMESTEPPING.MAXTIME",     po::value<double>(&maxTime), "Max Time")
	    ("TIMESTEPPING.MAXTIMESTEP", po::value<int>(&maxTimeStep), "Max Time Step")
	    ("TIMESTEPPING.FILTERSTEP",  po::value<int>(&filterStep), "Filter Step")
	    ("TIMESTEPPING.CHECKSTEP",   po::value<int>(&checkStep), "Check Step")
	    ("TIMESTEPPING.DUMPSTEP",    po::value<int>(&dumpStep), "Dump Step")
	    ("TIMESTEPPING.IMAGESTEP",   po::value<int>(&imageStep), "Image Step")
	    ("BC.BCXTYPE", 		 po::value<string>(&bcXType_str), "BC Type in X") 
	    ("BC.BCYTYPE", 		 po::value<string>(&bcYType_str), "BC Type in Y") 
	    ("BC.BCZTYPE", 		 po::value<string>(&bcZType_str), "BC Type in Z") 
	    ("BC.BCX0",			 po::value<string>(&bcX0_str), "BC Kind for X0")
	    ("BC.BCX1",			 po::value<string>(&bcX1_str), "BC Kind for X1")
	    ("BC.BCY0"	,		 po::value<string>(&bcY0_str), "BC Kind for Y0")
	    ("BC.BCY1",			 po::value<string>(&bcY1_str), "BC Kind for Y1")
	    ("BC.BCZ0",			 po::value<string>(&bcZ0_str), "BC Kind for Z0")
	    ("BC.BCZ1",			 po::value<string>(&bcZ1_str), "BC Kind for Z1")
	    ("DOMAIN.NX",                po::value<int>(&Nx), "Number of points in X-Direction")
	    ("DOMAIN.NY",                po::value<int>(&Ny), "Number of points in Y-Direction")
	    ("DOMAIN.NZ",                po::value<int>(&Nz), "Number of points in Z-Direction")
	    ("PENCILDECOMP.PROW",        po::value<int>(&pRow), "Number of Rows in the decomposition")
	    ("PENCILDECOMP.PCOL",        po::value<int>(&pCol), "Number of Cols in the decomposition")
	    ("SOLVER.ALPHAF", 		 po::value<double>(&alphaF), "Filter alpha coefficient value")
	    ("SOLVER.MU_REF", 		 po::value<double>(&mu_ref), "Reference Viscosity")
	    ("SOLVER.USETIMING",         po::value<bool>(&useTiming), "Report timing for different code sections")
	    ("RESTART.FROMRESTART", 	 po::value<bool>(&fromRestart), "Do we start the simulation from a restart")
	    ("RESTART.FILENAME",	 po::value<string>(&filename), "Filename of the restart file");
	
	po::store(po::parse_config_file(input_file, input), vm);    
	po::notify(vm);

	//Running through all of the input parameters
	
	cout << " ============== " << endl;
	cout << " =INPUT PARAMS= " << endl;
	cout << " ============== " << endl;

	//Should also probably check that these are positive

	//MISSING A THING FOR WHAT KIND OF TIMESTEPPING IT IS (CONST CFL OR DT)
	checkValue<double>("TIMESTEPPING.CFL", "CFL", CFL, 0.5);
	checkValue<double>("TIMESTEPPING.MAXTIME", "maxTime", maxTime, 1000);
	checkValue<int>(   "TIMESTEPPING.MAXTIMESTEP", "maxTimeStep", maxTimeStep, 10000);
	checkValue<int>(   "TIMESTEPPING.FILTERSTEP", "filterStep", filterStep, 1);
	checkValue<int>(   "TIMESTEPPING.CHECKSTEP", "checkStep", checkStep, 1);
	checkValue<int>(   "TIMESTEPPING.DUMPSTEP", "dumpStep", dumpStep, 1000);
	checkValue<int>(   "TIMESTEPPING.IMAGESTEP", "imageStep", imageStep, 25); //IS THIS OBSOLITE?

	forceValue<int>("DOMAIN.NX", "Nx", Nx);
	forceValue<int>("DOMAIN.NY", "Ny", Ny);
	forceValue<int>("DOMAIN.NZ", "Nz", Nz);

	checkValue<int>("PENCILDECOMP.PROW", "pRow", pRow, 0);
	checkValue<int>("PENCILDECOMP.PCOL", "pCol", pCol, 0);

	checkValue<double>("SOLVER.ALPHAF", "alphaF", alphaF, 0.45);
	forceValue<double>("SOLVER.MU_REF", "mu_ref", mu_ref);
	checkValue<bool>("SOLVER.USETIMING", "useTiming", useTiming, false);

	checkValue<bool>("RESTART.FROMRESTART", "fromRestart", fromRestart, false);
	if(fromRestart == true){
	    forceValue<string>("RESTART.FILENAME", "filename", filename);
	}
	
	//Parse from string...
	parseBCTypeFromString("BC.BCXTYPE", bcXType_str, bcXType);
	parseBCTypeFromString("BC.BCYTYPE", bcYType_str, bcYType);
	parseBCTypeFromString("BC.BCZTYPE", bcZType_str, bcZType);
	parseBCKindFromString("BC.BCX0", bcX0_str, bcX0);
	parseBCKindFromString("BC.BCX1", bcX1_str, bcX1);
	parseBCKindFromString("BC.BCY0", bcY0_str, bcY0);
	parseBCKindFromString("BC.BCY1", bcY1_str, bcY1);
	parseBCKindFromString("BC.BCZ0", bcZ0_str, bcZ0);
	parseBCKindFromString("BC.BCZ1", bcZ1_str, bcZ1);

	//NEED TO ADD IN PERIODIC DISPLACEMENTS...

	//FROM RESTART STUFF DOES NOTHING RIGHT NOW...
      }

      //Now we'll have to broadcast all of that stuff out to the other ranks...
      MPI_Bcast(&CFL, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

    }	

    template<typename T>
    void checkValue(string vmKey, string informalName, T inputContainer, T defaultVal){

	if(vm.count(vmKey)){
	    cout << " > " << informalName << " = " << inputContainer << endl;
	}else{
	    cout << " > " << informalName << " not specified, using default value of " << defaultVal << endl;
	    inputContainer = defaultVal;
	}
    }

    template<typename T>
    void forceValue(string vmKey, string informalName, T inputContainer){

	if(vm.count(vmKey)){
	    cout << " > " << informalName << " = " << inputContainer << endl;
	}else{
	    cout << " > " << informalName << " not specified, NEEDS TO BE SPECIFED " << endl;
	    MPI_Abort(MPI_COMM_WORLD, -10);	    
	}
    }

    void parseBCTypeFromString(string vmKey, string inString, BC::BCType currentType){

	if(vm.count(vmKey)){
	    
	    if(strcmp(inString.c_str(), "PERIODIC_SOLVE")==0){
		currentType = BC::PERIODIC_SOLVE;
	    }else if(strcmp(inString.c_str(), "DIRICHLET_SOLVE")==0){
	 	currentType = BC::DIRICHLET_SOLVE;
	    }else{
		cout << " > UNKNOWN BOUNDARY TYPE SPECIFIED: " << inString << endl;
		MPI_Abort(MPI_COMM_WORLD, -10);
	    }

	    cout << " > " << vmKey << " = " << inString << endl;
	    
	}else{
	    cout << " > " << vmKey << " = " << inString << " not specified " << endl;
	    MPI_Abort(MPI_COMM_WORLD, -10);	    
	}
    }

    void parseBCKindFromString(string vmKey, string inString, BC::BCKind currentType){

	if(vm.count(vmKey)){
	    
	    if(strcmp(inString.c_str(), "INTERNALLY_PERIODIC")==0){
		currentType = BC::INTERNALLY_PERIODIC;
	    }else if(strcmp(inString.c_str(), "PERIODIC")==0){
	 	currentType = BC::PERIODIC;
	    }else if(strcmp(inString.c_str(), "ADIABATIC_WALL")==0){
	 	currentType = BC::ADIABATIC_WALL;
	    }else if(strcmp(inString.c_str(), "SPONGE")==0){
	 	currentType = BC::SPONGE;
	    }else if(strcmp(inString.c_str(), "RECT_CURVILINEARSPONGE")==0){
	 	currentType = BC::RECT_CURVILINEARSPONGE;
	    }else if(strcmp(inString.c_str(), "CYL_CURVILINEARSPONGE")==0){
	 	currentType = BC::CYL_CURVILINEARSPONGE;
	    }else if(strcmp(inString.c_str(), "SPHERICAL_CURVILINEARSPONGE")==0){
	 	currentType = BC::SPHERICAL_CURVILINEARSPONGE;
	    }else if(strcmp(inString.c_str(), "CONST_T_WALL")==0){
	 	currentType = BC::CONST_T_WALL;
	    }else if(strcmp(inString.c_str(), "MOVING_ADIABATIC_WALL")==0){
	 	currentType = BC::MOVING_ADIABATIC_WALL;
	    }else if(strcmp(inString.c_str(), "INLET")==0){
	 	currentType = BC::INLET;
	    }else{
		cout << " > UNKNOWN BOUNDARY TYPE SPECIFIED: " << inString << endl;
		MPI_Abort(MPI_COMM_WORLD, -10);
	    }

	    cout << " > " << vmKey << " = " << inString << endl;

	}else{
	    cout << " > " << vmKey << " = " << inString << " not specified" << endl;
	    MPI_Abort(MPI_COMM_WORLD, -10);	    
	}
    }


};

#endif 
