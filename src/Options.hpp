#ifndef _OPTIONSH_
#define _OPTIONSH_

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <boost/program_options.hpp>

#include "TimeStepping.hpp"
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
	TimeStepping::TimeSteppingType TSType;
	string TSType_str;
	double CFL, dt, maxTime;
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
	    ("TIMESTEPPING.TSTYPE",	 po::value<string>(&TSType_str), "Time Stepping Type")
	    ("TIMESTEPPING.CFL", 	 po::value<double>(&CFL), "CFL Number")
	    ("TIMESTEPPING.DT", 	 po::value<double>(&dt), "Time Step Size")
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
	    ("BC.PERIODICDISPX_X", 	 po::value<double>(&periodicDisp[0][0]), "Periodic Displacement for X0 face to X1 face in the x-direction")
	    ("BC.PERIODICDISPX_Y", 	 po::value<double>(&periodicDisp[0][1]), "Periodic Displacement for X0 face to X1 face in the y-direction")
	    ("BC.PERIODICDISPX_Z", 	 po::value<double>(&periodicDisp[0][2]), "Periodic Displacement for X0 face to X1 face in the z-direction")
	    ("BC.PERIODICDISPY_X", 	 po::value<double>(&periodicDisp[1][0]), "Periodic Displacement for Y0 face to Y1 face in the x-direction")
	    ("BC.PERIODICDISPY_Y", 	 po::value<double>(&periodicDisp[1][1]), "Periodic Displacement for Y0 face to Y1 face in the y-direction")
	    ("BC.PERIODICDISPY_Z", 	 po::value<double>(&periodicDisp[1][2]), "Periodic Displacement for Y0 face to Y1 face in the z-direction")
	    ("BC.PERIODICDISPZ_X", 	 po::value<double>(&periodicDisp[2][0]), "Periodic Displacement for Z0 face to Z1 face in the x-direction")
	    ("BC.PERIODICDISPZ_Y", 	 po::value<double>(&periodicDisp[2][1]), "Periodic Displacement for Z0 face to Z1 face in the y-direction")
	    ("BC.PERIODICDISPZ_Z", 	 po::value<double>(&periodicDisp[2][2]), "Periodic Displacement for Z0 face to Z1 face in the z-direction")
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
	parseTSTypeFromString("TIMESTEPPING.TSTYPE", TSType_str, TSType);
	if(TSType == TimeStepping::CONST_DT){
	    forceValue<double>("TIMESTEPPING.DT", "DT", dt);
	    CFL = -1.0;
	}else if(TSType == TimeStepping::CONST_CFL){
	    checkValue<double>("TIMESTEPPING.CFL", "CFL", CFL, 0.5);
	    dt = -1.0;
	}

	checkValue<double>("TIMESTEPPING.MAXTIME", "maxTime", maxTime, 1000);
	checkValue<int>(   "TIMESTEPPING.MAXTIMESTEP", "maxTimeStep", maxTimeStep, 10000);
	checkValue<int>(   "TIMESTEPPING.FILTERSTEP", "filterStep", filterStep, 1);
	checkValue<int>(   "TIMESTEPPING.CHECKSTEP", "checkStep", checkStep, 1);
	checkValue<int>(   "TIMESTEPPING.DUMPSTEP", "dumpStep", dumpStep, 1000);
	checkValue<int>(   "TIMESTEPPING.IMAGESTEP", "imageStep", imageStep, 25); //IS THIS OBSOLITE?

	forceValue<int>("DOMAIN.NX", "Nx", Nx);
	forceValue<int>("DOMAIN.NY", "Ny", Ny);
	forceValue<int>("DOMAIN.NZ", "Nz", Nz);

	//As of right now we're restricted to Nx|Ny|Nz > 10...
	if(Nx <= 10 || Ny <= 10 || Nz <= 10){
	    cout << " > Number of grid points in any direction must be greater than 10! " << endl;
	    MPI_Abort(MPI_COMM_WORLD, -10);	    
	}

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

	//PERIODIC BC DISPLACEMENTS...
	if(bcX0 == BC::PERIODIC && bcX1 == BC::PERIODIC){
	    forceValue<double>("BC.PERIODICDISPX_X", "x0 to x1 periodic displacement in x", periodicDisp[0][0]); 
	    forceValue<double>("BC.PERIODICDISPX_Y", "x0 to x1 periodic displacement in y", periodicDisp[0][1]); 
	    forceValue<double>("BC.PERIODICDISPX_Z", "x0 to x1 periodic displacement in z", periodicDisp[0][2]); 
	}else if((bcX0 == BC::PERIODIC && bcX1 != BC::PERIODIC) || (bcX0 != BC::PERIODIC && bcX1 == BC::PERIODIC)){
	    cout << " > " << " BC Mismatch for X0-X1 Faces! " << endl;
	    MPI_Abort(MPI_COMM_WORLD, -10);	    
	}  

	if(bcY0 == BC::PERIODIC && bcY1 == BC::PERIODIC){
	    forceValue<double>("BC.PERIODICDISPY_X", "y0 to y1 periodic displacement in x", periodicDisp[1][0]); 
	    forceValue<double>("BC.PERIODICDISPY_Y", "y0 to y1 periodic displacement in y", periodicDisp[1][1]); 
	    forceValue<double>("BC.PERIODICDISPY_Z", "y0 to y1 periodic displacement in z", periodicDisp[1][2]); 
	}else if((bcY0 == BC::PERIODIC && bcY1 != BC::PERIODIC) || (bcY0 != BC::PERIODIC && bcY1 == BC::PERIODIC)){
	    cout << " > " << " BC Mismatch for Y0-Y1 Faces! " << endl;
	    MPI_Abort(MPI_COMM_WORLD, -10);	    
	}  

	if(bcZ0 == BC::PERIODIC && bcZ1 == BC::PERIODIC){
	    forceValue<double>("BC.PERIODICDISPZ_X", "z0 to z1 periodic displacement in x", periodicDisp[2][0]); 
	    forceValue<double>("BC.PERIODICDISPZ_Y", "z0 to z1 periodic displacement in y", periodicDisp[2][1]); 
	    forceValue<double>("BC.PERIODICDISPZ_Z", "z0 to z1 periodic displacement in z", periodicDisp[2][2]); 
	}else if((bcZ0 == BC::PERIODIC && bcZ1 != BC::PERIODIC) || (bcZ0 != BC::PERIODIC && bcZ1 == BC::PERIODIC)){
	    cout << " > " << " BC Mismatch for Z0-Z1 Faces! " << endl;
	    MPI_Abort(MPI_COMM_WORLD, -10);	    
	}  

	//Do boundary condition validation to make sure things are peachy...

	//FROM RESTART STUFF DOES NOTHING RIGHT NOW...
      }

      //Now we'll have to broadcast all of that stuff out to the other ranks...
      MPI_Bcast(&CFL, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

    }	

    template<typename T>
    void checkValue(string vmKey, string informalName, T &inputContainer, T defaultVal){

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

    void parseTSTypeFromString(string vmKey, string inString, TimeStepping::TimeSteppingType &currentType);
    void parseBCTypeFromString(string vmKey, string inString, BC::BCType &currentType);
    void parseBCKindFromString(string vmKey, string inString, BC::BCKind &currentType);
    void bcValidation();

};

#endif 
