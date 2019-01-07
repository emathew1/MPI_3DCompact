#ifndef _OPTIONSH_
#define _OPTIONSH_

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>

#include "BC.hpp"

using namespace std;
namespace po = boost::program_options;

class Options{

    public:
	
	//Variable map...
	po::variable_map vm;

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
    
    Options(){
	
	ifstream input_file("MPI_3DCompact.in");

	po::options_description input("Input Parameters");

	input.add_options()
	    ("TIMESTEPPING.CFL", 	 po::value<double(&CFL), "CFL Number")
	    ("TIMESTEPPING.MAXTIME",     po::value<double(&maxTime), "Max Time")
	    ("TIMESTEPPING.MAXTIMESTEP", po::value<int>(&maxTimeStep), "Max Time Step")
	    ("TIMESTEPPING.FILTERSTEP",  po::value<int>(&filterStep), "Filter Step")
	    ("TIMESTEPPING.CHECKSTEP",   po::value<int>(&checkStep), "Check Step")
	    ("TIMESTEPPING.DUMPSTEP",    po::value<int>(&dumpStep), "Dump Step")
	    ("TIMESTEPPING.IMAGESTEP",   po::value<int>(&imageStep), "Image Step")
	    ("BC.BCXTYPE", 		 po::value<string>(&bcXType_str), "BC Type in X") 
	    ("BC.BCYTYPE", 		 po::value<string>(&bcYType_str), "BC Type in Y") 
	    ("BC.BCZTYPE", 		 po::value<string>(&bcZType_str), "BC Type in Z") 
	    ("BC.BCX0",			 po::value<string>(&bX0_str), "BC Kind for X0")
	    ("BC.BCX1",			 po::value<string>(&bX1_str), "BC Kind for X1")
	    ("BC.BCY0"	,		 po::value<string>(&bY0_str), "BC Kind for Y0")
	    ("BC.BCY1",			 po::value<string>(&bY1_str), "BC Kind for Y1")
	    ("BC.BCZ0",			 po::value<string>(&bZ0_str), "BC Kind for Z0")
	    ("BC.BCZ1",			 po::value<string>(&bZ1_str), "BC Kind for Z1")
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

	checkValue<double>("TIMESTEPPING.CFL", "CFL", CFL, 0.5);
	checkValue<double>("TIMESTEPPING.MAXTIME", "maxTime", maxTime, 1000);
	checkValue<int>(   "TIMESTEPPING.MAXTIMESTEP", "maxTimeStep", maxTimeStep, 10000);
	checkValue<int>(   "TIMESTEPPING.FILTERSTEP", "filterStep", filterStep, 1);
	checkValue<int>(   "TIMESTEPPING.CHECKSTEP", "checkStep", checkStep, 1);
	checkValue<int>(   "TIMESTEPPING.DUMPSTEP", "dumpStep", dumpStep, 1000);
	checkValue<int>(   "TIMESTEPPING.IMAGESTEP", "imageStep", imageStep, 25);

	forceValue<int>("DOMAIN.NX", "Nx");
	forceValue<int>("DOMAIN.NY", "Ny");
	forceValue<int>("DOMAIN.NZ", "Nz");

	checkValue<int>("PENCILDECOMP.PROW", "pRow", pRow, 0);
	checkValue<int>("PENCILDECOMP.PCOL", "pCol", pCol, 0);

	checkValue<double>("SOLVER.ALPHAF", "alphaF", alphaF, 0.45);
	forceValue<double>("SOLVER.MU_REF", "mu_ref");
	checkValue<bool>("SOLVER.USETIME", "useTiming", useTiming, false);

	checkValue<bool>("RESTART.FROMRESTART", "fromRestart", fromRestart, false);
	if(fromRestart == true){
	    forceValue<string>("RESTART.FILENAME", "filename");
	}
	
	//Parse from string...
	parseBCTypeFromString("BC.BCXTYPE", bcXType_str, bcXType);
	parseBCTypeFromString("BC.BCYTYPE", bcYType_str, bcYType);
	parseBCTypeFromString("BC.BCZTYPE", bcZType_str, bcZType);

	//Need to parse the other boundary stuff...

    }	

    template <T>;
    void checkValue(string vmKey, string informalName, T inputContainer, T defaultVal){

	if(vm.count(vmKey)){
	    cout << " > " << informalName << " = " << inputContainer << endl;
	}else{
	    cout << " > " << informalName << " not specified, using default value of " << defaultVal << endl;
	    inputContainer = defaultVal;
	}
    }

    void forceValue(string vmKey, string informalName){

	if(vm.count(vmKey)){
	    cout << " > " << informalName << " = " << inputContainer << endl;
	}else{
	    cout << " > " << informalName << " not specified, NEEDS TO BE SPECIFED " << defaultVal << endl;
	    MPI_Abort(MPI_COMM_WORLD, -10);	    
	}
    }

    void parseBCTypeFromString(string vmKey, string inString, BCType currentType){

	if(vm.count(vmKey){
	    
	    if(strcmp(inString, "PERIODIC_SOLVE")==0){
		currentType = BC::PERIODIC_SOLVE;
	    }else if(strcmp(inString, "DIRICHLET_SOLVE")==0){
	 	currentType = BC::DIRICHLET_SOLVE;
	    }else{
		cout << " > UNKNOWN BOUNDARY TYPE SPECIFIED: " << inString << endl;
		MPI_Abort(MPI_COMM_WORLD, -10);
	    }
	}else{
	    cout << " > " << informalName << " not specified, NEEDS TO BE SPECIFED " << defaultVal << endl;
	    MPI_Abort(MPI_COMM_WORLD, -10);	    
	}
    }

};

#endif 
