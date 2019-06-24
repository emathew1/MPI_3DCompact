#ifndef _OPTIONSH_
#define _OPTIONSH_

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

class Options{

    public:

	//Moving all enums here so Option.hpp can be non-dependent on anything else...
        enum TimeSteppingType {CONST_DT, CONST_CFL};
        enum BCType {PERIODIC_SOLVE, DIRICHLET_SOLVE};
        enum BCKind {PERIODIC,
                     INTERNALLY_PERIODIC,
                     SPONGE,
                     ADIABATIC_WALL,
                     CONST_T_WALL,
                     MOVING_ADIABATIC_WALL,
                     INLET};
        enum SpongeKind {RECTILINEAR, CYLINDRICAL};
	
	enum FDType {CD2, PADE6, PENTA10};
	enum FilterType {COMPACT8, COMPACT10};

	enum LESModel {NOMODEL, VREMAN, DSM};

	enum StatsAvgType {NONE, XI1_AVG, XI2_AVG, XI3_AVG, LOCAL, GLOBAL};

	enum RKType {TVDRK3, RK4, KENRK4, LSLDDRK4};

	//Variable map...
	po::variables_map vm;

	//Root rank
	int root;
	int mpiRank;

	//Time Stepping Information
	TimeSteppingType TSType;
	string TSType_str;
	double CFL, dt, maxTime;
	int maxTimeStep, filterStep, checkStep, dumpStep;
	bool subStepFiltering; 

	//Boundary Condition Stuff
	string bcXType_str, bcYType_str, bcZType_str;
	BCType bcXType, bcYType, bcZType;

	string bcX0_str, bcX1_str, bcY0_str, bcY1_str, bcZ0_str, bcZ1_str;
	BCKind bcX0, bcX1, bcY0, bcY1, bcZ0, bcZ1;
	double periodicDisp[3][3];

	//SPONGE STUFF
	string spongeKind_str;
	SpongeKind spongeKind;
	double spongeP, spongeAvgT, spongeStrength;
	double spongeRectX0Perc, spongeRectY0Perc, spongeRectZ0Perc;
	double spongeRectX1Perc, spongeRectY1Perc, spongeRectZ1Perc;
	int spongeCylAxisOrient;
	double spongeCylAxisX, spongeCylAxisY, spongeCylAxisZ, rMin;  

	//DOMAIN Information
	int Nx, Ny, Nz;

	//Pencil Decomposition Stuff
	int pRow, pCol;

	//SOLVER Information
	double alphaF, mu_ref;
	bool useTiming;
	FDType xFDType, yFDType, zFDType;
	string xFDType_str, yFDType_str, zFDType_str;
	FilterType filterType;
	string filterType_str;
	RKType rkType;
	string rkType_str;	

	//RESTART stuff
	bool fromRestart;
	bool onlyGridFromRestart;
	string filename;
	bool spongeFromRestart;
	string sponge_filename;
	int  timeStep;    
	double time;


	//STATS stuff
	StatsAvgType statsAvgType;
	string statsAvgType_str;

	bool velocityStats;
	bool velocityStatsFromRestart;
	string velocityStatsFilename;

	bool thermoStats;
	bool thermoStatsFromRestart;
	string thermoStatsFilename;

	int stats_interval;

	//LES stuff
	StatsAvgType musgsAvgType;
	string musgsAvgType_str;
	LESModel lesModel;
	string lesModel_str;	

	//DSM LES Model Stuff
	bool dumpCoeffRange;
	bool useTaukk;
	int testFilterType;

    Options(int mpiRank){
   

      root = 0;
      this->mpiRank = mpiRank;

 
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
	    ("TIMESTEPPING.SUBSTEPFILTERING",  po::value<bool>(&subStepFiltering), "Filtering every substep")
	    ("TIMESTEPPING.CHECKSTEP",   po::value<int>(&checkStep), "Check Step")
	    ("TIMESTEPPING.DUMPSTEP",    po::value<int>(&dumpStep), "Dump Step")
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
	    ("SPONGE.KIND",		 po::value<string>(&spongeKind_str), "Type of sponge boundary")
	    ("SPONGE.AVGTIME",		 po::value<double>(&spongeAvgT), "Sponge averaging time")
	    ("SPONGE.PINF",		 po::value<double>(&spongeP), "Sponge back pressure")
	    ("SPONGE.STRENGTH",		 po::value<double>(&spongeStrength), "Sponge strength")
	    ("SPONGE.RECTX0PERC",	 po::value<double>(&spongeRectX0Perc), "Sponge Percent of Domain in X Direction - X0 Face")
	    ("SPONGE.RECTY0PERC",	 po::value<double>(&spongeRectY0Perc), "Sponge Percent of Domain in Y Direction - Y0 Face")
	    ("SPONGE.RECTZ0PERC",	 po::value<double>(&spongeRectZ0Perc), "Sponge Percent of Domain in Z Direction - Z0 Face")
	    ("SPONGE.RECTX1PERC",	 po::value<double>(&spongeRectX1Perc), "Sponge Percent of Domain in X Direction - X1 Face")
	    ("SPONGE.RECTY1PERC",	 po::value<double>(&spongeRectY1Perc), "Sponge Percent of Domain in Y Direction - Y1 Face")
	    ("SPONGE.RECTZ1PERC",	 po::value<double>(&spongeRectZ1Perc), "Sponge Percent of Domain in Z Direction - Z1 Face")
    	    ("SPONGE.CYLAXIS_ORIENT",	 po::value<int>(&spongeCylAxisOrient), "Axis that the cylindrical sponge shape faces")
	    ("SPONGE.CYLAXIS_X",	 po::value<double>(&spongeCylAxisX), "X-location of cylindrical sponge center")
	    ("SPONGE.CYLAXIS_Y", 	 po::value<double>(&spongeCylAxisY), "Y-location of cylindrical sponge center")
	    ("SPONGE.CYLAXIS_Z",	 po::value<double>(&spongeCylAxisZ), "Z-location of cylindrical sponge center")	
	    ("SPONGE.RMIN",		 po::value<double>(&rMin), "Minimum radius of cylinder sponge")	
	    ("DOMAIN.NX",                po::value<int>(&Nx), "Number of points in X-Direction")
	    ("DOMAIN.NY",                po::value<int>(&Ny), "Number of points in Y-Direction")
	    ("DOMAIN.NZ",                po::value<int>(&Nz), "Number of points in Z-Direction")
	    ("PENCILDECOMP.PROW",        po::value<int>(&pRow), "Number of Rows in the decomposition")
	    ("PENCILDECOMP.PCOL",        po::value<int>(&pCol), "Number of Cols in the decomposition")
	    ("SOLVER.ALPHAF", 		 po::value<double>(&alphaF), "Filter alpha coefficient value")
	    ("SOLVER.MU_REF", 		 po::value<double>(&mu_ref), "Reference Viscosity")
	    ("SOLVER.USETIMING",         po::value<bool>(&useTiming), "Report timing for different code sections")
	    ("SOLVER.XFDTYPE",		 po::value<string>(&xFDType_str), "Finite Differences Type in the X-Direction")
	    ("SOLVER.YFDTYPE",		 po::value<string>(&yFDType_str), "Finite Differences Type in the Y-Direction")
	    ("SOLVER.ZFDTYPE",		 po::value<string>(&zFDType_str), "Finite Differences Type in the Z-Direction")
	    ("SOLVER.FILTERTYPE",	 po::value<string>(&filterType_str), "Filter type for all three directions")
	    ("SOLVER.RKTYPE",		 po::value<string>(&rkType_str), "Type of Runge-Kutta time stepping")
	    ("RESTART.FROMRESTART", 	 po::value<bool>(&fromRestart), "Do we start the simulation from a restart")
	    ("RESTART.ONLYGRIDFROMRESTART", 	 po::value<bool>(&onlyGridFromRestart), "Only pull the grid from the restart file")
	    ("RESTART.FILENAME",	 po::value<string>(&filename), "Filename of the restart file")
	    ("RESTART.SPONGEAVGFROMRESTART", 	 po::value<bool>(&spongeFromRestart), "Pull the sponge average from file")
	    ("RESTART.SPONGEAVGFILENAME",	 po::value<string>(&sponge_filename), "Filename of the sponge-avg file")
	    ("STATS.VELOCITYSTATS",	 po::value<bool>(&velocityStats), "Flag for recording velocity statistics")
	    ("STATS.THERMOSTATS",	 po::value<bool>(&thermoStats), "Flag for recording thermodynamic quantity statistics")
	    ("STATS.VELOCITYSTATSFROMRESTART",	 po::value<bool>(&velocityStatsFromRestart), "Flag for using restart in velocity stats")
	    ("STATS.THERMOSTATSFROMRESTART",	 po::value<bool>(&thermoStatsFromRestart), "Flag for using restart in thermo statistics")
	    ("STATS.VELOCITYSTATSFILENAME",	 po::value<string>(&velocityStatsFilename), "Velocity statistics restart filename")
	    ("STATS.THERMOSTATSFILENAME",	 po::value<string>(&thermoStatsFilename), "Thermo statistics restart filename")
	    ("STATS.STATSAVGTYPE",	 po::value<string>(&statsAvgType_str), "Name of spacial averaging type for statistics")
	    ("STATS.INTERVAL",	 po::value<int>(&stats_interval), "How often (in terms of timesteps) do we update the statistics data")
	    ("LES.LESMODEL",	 po::value<string>(&lesModel_str), "Name of the LES model")
	    ("LES.LESAVERAGING", po::value<string>(&musgsAvgType_str), "Averaging of mu_sgs type")
	    ("LES.DSMTESTFILTER", po::value<int>(&testFilterType), "Type of test filter used in dynamic coefficient calculation")
	    ("LES.DSMUSETAUKK", po::value<bool>(&useTaukk), "Whether or not to use taukk component in pressure calculation")
	    ("LES.DSMCOEFFRANGE", po::value<bool>(&dumpCoeffRange), "Whether or not to dump the Smag. Coeff, CI, and Prt.");

	
	    //Potentially include the style of computation options from CurvilinearCsolver? OCC vs. VANILLA etc.	


	po::store(po::parse_config_file(input_file, input), vm);    
	po::notify(vm);

	//Running through all of the input parameters
	
	cout << " ============== " << endl;
	cout << " =INPUT PARAMS= " << endl;
	cout << " ============== " << endl;

	//Should also probably check that these are positive

	//MISSING A THING FOR WHAT KIND OF TIMESTEPPING IT IS (CONST CFL OR DT)
	parseTSTypeFromString("TIMESTEPPING.TSTYPE", TSType_str, TSType);
	if(TSType == CONST_DT){
	    forceValue<double>("TIMESTEPPING.DT", "DT", dt);
	    CFL = -1.0;
	}else if(TSType == CONST_CFL){
	    checkValue<double>("TIMESTEPPING.CFL", "CFL", CFL, 0.5);
	    dt = -1.0;
	}

	checkValue<double>("TIMESTEPPING.MAXTIME", "maxTime", maxTime, 1000);
	checkValue<int>(   "TIMESTEPPING.MAXTIMESTEP", "maxTimeStep", maxTimeStep, 10000);
	checkValue<int>(   "TIMESTEPPING.FILTERSTEP", "filterStep", filterStep, 1);
	checkValue<int>(   "TIMESTEPPING.CHECKSTEP", "checkStep", checkStep, 1);
	checkValue<int>(   "TIMESTEPPING.DUMPSTEP", "dumpStep", dumpStep, 1000);
	checkValue<bool>(   "TIMESTEPPING.SUBSTEPFILTERING", "subStepFiltering", subStepFiltering, false);

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
	parseFDTypeFromString("SOLVER.XFDTYPE", xFDType_str, xFDType);
	parseFDTypeFromString("SOLVER.YFDTYPE", yFDType_str, yFDType);
	parseFDTypeFromString("SOLVER.ZFDTYPE", zFDType_str, zFDType);
	parseFilterTypeFromString("SOLVER.FILTERTYPE", filterType_str, filterType);
	parseRKTypeFromString("SOLVER.RKTYPE", rkType_str, rkType);

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
	if(bcX0 == PERIODIC && bcX1 == PERIODIC){
	    forceValue<double>("BC.PERIODICDISPX_X", "x0 to x1 periodic displacement in x", periodicDisp[0][0]); 
	    forceValue<double>("BC.PERIODICDISPX_Y", "x0 to x1 periodic displacement in y", periodicDisp[0][1]); 
	    forceValue<double>("BC.PERIODICDISPX_Z", "x0 to x1 periodic displacement in z", periodicDisp[0][2]); 
	}else if((bcX0 == PERIODIC && bcX1 != PERIODIC) || (bcX0 != PERIODIC && bcX1 == PERIODIC)){
	    cout << " > " << " BC Mismatch for X0-X1 Faces! " << endl;
	    MPI_Abort(MPI_COMM_WORLD, -10);	    
	}else{
	    periodicDisp[0][0] = 0.0;
	    periodicDisp[0][1] = 0.0;
	    periodicDisp[0][2] = 0.0;
	}  

	if(bcY0 == PERIODIC && bcY1 == PERIODIC){
	    forceValue<double>("BC.PERIODICDISPY_X", "y0 to y1 periodic displacement in x", periodicDisp[1][0]); 
	    forceValue<double>("BC.PERIODICDISPY_Y", "y0 to y1 periodic displacement in y", periodicDisp[1][1]); 
	    forceValue<double>("BC.PERIODICDISPY_Z", "y0 to y1 periodic displacement in z", periodicDisp[1][2]); 
	}else if((bcY0 == PERIODIC && bcY1 != PERIODIC) || (bcY0 != PERIODIC && bcY1 == PERIODIC)){
	    cout << " > " << " BC Mismatch for Y0-Y1 Faces! " << endl;
	    MPI_Abort(MPI_COMM_WORLD, -10);	    
	}else{
	    periodicDisp[1][0] = 0.0;
	    periodicDisp[1][1] = 0.0;
	    periodicDisp[1][2] = 0.0;
	}  

	if(bcZ0 == PERIODIC && bcZ1 == PERIODIC){
	    forceValue<double>("BC.PERIODICDISPZ_X", "z0 to z1 periodic displacement in x", periodicDisp[2][0]); 
	    forceValue<double>("BC.PERIODICDISPZ_Y", "z0 to z1 periodic displacement in y", periodicDisp[2][1]); 
	    forceValue<double>("BC.PERIODICDISPZ_Z", "z0 to z1 periodic displacement in z", periodicDisp[2][2]); 
	}else if((bcZ0 == PERIODIC && bcZ1 != PERIODIC) || (bcZ0 != PERIODIC && bcZ1 == PERIODIC)){
	    cout << " > " << " BC Mismatch for Z0-Z1 Faces! " << endl;
	    MPI_Abort(MPI_COMM_WORLD, -10);	    
	}else{
	    periodicDisp[2][0] = 0.0;
	    periodicDisp[2][1] = 0.0;
	    periodicDisp[2][2] = 0.0;
	}  

	//Do boundary condition validation to make sure things are peachy...
	bcValidation();

	//Lets collect and parse all of the sponge stuff...
	if(bcX0 == SPONGE || \
	   bcX1 == SPONGE || \
	   bcY0 == SPONGE || \
	   bcY1 == SPONGE || \
	   bcZ0 == SPONGE || \
	   bcZ1 == SPONGE){

	    parseSpongeFromString("SPONGE.KIND", spongeKind_str, spongeKind); 
	    checkValue<double>("SPONGE.AVGTIME",  "spongeAvgT", spongeAvgT, 1.0);
	    checkValue<double>("SPONGE.PINF",     "spongeP"   , spongeP, 1.0/1.4);
	    checkValue<double>("SPONGE.STRENGTH", "spongeStrength", spongeStrength, 1.0);

	    if(spongeKind == RECTILINEAR){

		if(bcX0 == SPONGE){
	    	    forceValue<double>("SPONGE.RECTX0PERC", "spongeRectX0Perc", spongeRectX0Perc); 
		}	

		if(bcX1 == SPONGE){
	    	    forceValue<double>("SPONGE.RECTX1PERC", "spongeRectX1Perc", spongeRectX1Perc); 
		}	

		if(bcY0 == SPONGE){
	    	    forceValue<double>("SPONGE.RECTY0PERC", "spongeRectY0Perc", spongeRectY0Perc); 
		}	

		if(bcY1 == SPONGE){
	    	    forceValue<double>("SPONGE.RECTY1PERC", "spongeRectY1Perc", spongeRectY1Perc); 
		}	

		if(bcZ0 == SPONGE){
	    	    forceValue<double>("SPONGE.RECTZ0PERC", "spongeRectZ0Perc", spongeRectZ0Perc); 
		}	

		if(bcZ1 == SPONGE){
	    	    forceValue<double>("SPONGE.RECTZ1PERC", "spongeRectZ1Perc", spongeRectZ1Perc); 
		}	

	    }else if(spongeKind == CYLINDRICAL){
		forceValue<int>("SPONGE.CYLAXIS_ORIENT", "spongeCylAxisOrient", spongeCylAxisOrient);
                if(spongeCylAxisOrient < 0 || spongeCylAxisOrient > 2){
                    cout << " > Sponge cylinder orientation invalid number, needs to be 0, 1, 2. value = " << spongeCylAxisOrient << endl;
                    MPI_Abort(MPI_COMM_WORLD, -10);
                }

		forceValue<double>("SPONGE.CYLAXIS_X", "spongeCylAxisX", spongeCylAxisX);
		forceValue<double>("SPONGE.CYLAXIS_Y", "spongeCylAxisY", spongeCylAxisY);
		forceValue<double>("SPONGE.CYLAXIS_Z", "spongeCylAxisZ", spongeCylAxisZ);
		forceValue<double>("SPONGE.RMIN", "rMin", rMin);
	    }

	}


	//Restart file stuff
	forceValue<bool>("RESTART.FROMRESTART", "fromRestart", fromRestart);
	forceValue<bool>("RESTART.ONLYGRIDFROMRESTART", "onlyGridFromRestart", onlyGridFromRestart);
	
	if(fromRestart || onlyGridFromRestart){
	    forceValue<string>("RESTART.FILENAME", "restart filename", filename);
	}

	if(fromRestart && onlyGridFromRestart){
	    cout << "Both 'only use grid from restart file' and 'from restart file' checked off, these are conflicting!" << endl;
	    MPI_Abort(MPI_COMM_WORLD, -10);
	}

	//Need the dimensions of the grid if not from restart frile
	if(!(fromRestart || onlyGridFromRestart)){
	    forceValue<int>("DOMAIN.NX", "Nx", Nx);
	    forceValue<int>("DOMAIN.NY", "Ny", Ny);
	    forceValue<int>("DOMAIN.NZ", "Nz", Nz);
	}

	forceValue<bool>("RESTART.SPONGEAVGFROMRESTART", "spongeFromRestart", spongeFromRestart);
	if(spongeFromRestart){
	    forceValue<string>("RESTART.SPONGEAVGFILENAME", "sponge_filename", sponge_filename);
	}
	

	//Now lets do all of the statistics stuff	
	checkValue<bool>("STATS.VELOCITYSTATS", "velocityStats", velocityStats, false);
	
	if(velocityStats){
	    forceValue<bool>("STATS.VELOCITYSTATSFROMRESTART", "velocityStatsFromRestart", velocityStatsFromRestart);
	    if(velocityStatsFromRestart){
		forceValue<string>("STATS.VELOCITYSTATSFILENAME", "velocityStatsFilename", velocityStatsFilename);
	    }
	}

	checkValue<bool>("STATS.THERMOSTATS",   "thermoStats",   thermoStats,   false);
	if(thermoStats){
	    forceValue<bool>("STATS.THERMOSTATSFROMRESTART", "thermoStatsFromRestart", thermoStatsFromRestart);
	    if(thermoStatsFromRestart){
		forceValue<string>("STATS.THERMOSTATSFILENAME", "thermoStatsFilename", thermoStatsFilename);
	    }

	}

	//Do the averaging type check here
	//parseStatsAvgType()...
	parseStatsAvgTypeString("STATS.STATSAVGTYPE", statsAvgType_str, statsAvgType);

	//Is the restart redundant?	
	checkValue<int>("STATS.INTERVAL", "stats_interval", stats_interval, 5);

	
	//LES model stuff
	parseLESModelString("LES.LESMODEL", lesModel_str, lesModel);
	parseStatsAvgTypeString("LES.LESAVERAGING", musgsAvgType_str, musgsAvgType);

	checkValue<int>("LES.DSMTESTFILTER", "testFilterType", testFilterType, 1);
	checkValue<bool>("LES.DSMUSETAUKK", "useTaukk", useTaukk, true);
	checkValue<bool>("LES.DSMCOEFFRANGE", "dumpCoeffRange", dumpCoeffRange, false);

      }


      //Now we'll have to broadcast all of that stuff out to the other ranks...

      //Probably not the most efficient way of doing this but definitely the laziest

      //[TIMESTEPPING]
      MPI_Bcast(&TSType, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&CFL, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&dt,  1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&maxTime, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&maxTimeStep, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&filterStep, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&checkStep, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&dumpStep, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&subStepFiltering, 1, MPI_INT, root, MPI_COMM_WORLD);

      //[BC]
      MPI_Bcast(&bcXType, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&bcYType, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&bcZType, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&bcX0, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&bcX1, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&bcY0, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&bcY1, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&bcZ0, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&bcZ1, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(periodicDisp, 9, MPI_DOUBLE, root, MPI_COMM_WORLD);

      //[SPONGE]
      MPI_Bcast(&spongeKind, 1, MPI_INT, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeP, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeAvgT, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeStrength, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeRectX0Perc, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeRectX1Perc, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeRectY0Perc, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeRectY1Perc, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeRectZ0Perc, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeRectZ1Perc, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeCylAxisOrient, 1, MPI_INT, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeCylAxisX, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeCylAxisY, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeCylAxisZ, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 
      MPI_Bcast(&rMin, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 

      //[DOMAIN]
      MPI_Bcast(&Nx, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&Ny, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&Nz, 1, MPI_INT, root, MPI_COMM_WORLD);

      //[PENCILDECOMP]
      MPI_Bcast(&pRow, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&pCol, 1, MPI_INT, root, MPI_COMM_WORLD);

      //[SOLVER]
      MPI_Bcast(&alphaF, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&mu_ref, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&useTiming, 1, MPI_C_BOOL, root, MPI_COMM_WORLD);
      MPI_Bcast(&xFDType, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&yFDType, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&zFDType, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&filterType, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&rkType, 1, MPI_INT, root, MPI_COMM_WORLD);
   

      //[RESTART]
      MPI_Bcast(&fromRestart, 1, MPI_C_BOOL, root, MPI_COMM_WORLD); 
      MPI_Bcast(&onlyGridFromRestart, 1, MPI_C_BOOL, root, MPI_COMM_WORLD); 
      MPI_Bcast(&spongeFromRestart, 1, MPI_C_BOOL, root, MPI_COMM_WORLD);     

      //Restart file name 
      int stringSize; 
      if(mpiRank == root){
          stringSize = filename.size()+1;
      }
      MPI_Bcast(&stringSize, 1, MPI_INT, root, MPI_COMM_WORLD);

      char filename_c[stringSize];
      if(mpiRank == root){
          strcpy(filename_c,filename.c_str());
      }
 
      MPI_Bcast(filename_c, stringSize, MPI_CHAR, root, MPI_COMM_WORLD);
      filename.assign(filename_c, stringSize);

      //Sponge Average Restart File name
      if(mpiRank == root){
          stringSize = sponge_filename.size()+1;
      }
      MPI_Bcast(&stringSize, 1, MPI_INT, root, MPI_COMM_WORLD);

      char sfilename_c[stringSize];
      if(mpiRank == root){
          strcpy(sfilename_c,sponge_filename.c_str());
      }
 
      MPI_Bcast(sfilename_c, stringSize, MPI_CHAR, root, MPI_COMM_WORLD);
      sponge_filename.assign(sfilename_c, stringSize);


      //STATS Stuff


      MPI_Bcast(&velocityStats, 1, MPI_C_BOOL, root, MPI_COMM_WORLD);
      if(velocityStats){

	  //if we're doing velocity stats, see if we need to pull from restartt
	  MPI_Bcast(&velocityStatsFromRestart, 1, MPI_C_BOOL, root, MPI_COMM_WORLD);

	  //If we're pulling from restart, get the filename
	  if(velocityStatsFromRestart){

	      if(mpiRank == root){
          	  stringSize = velocityStatsFilename.size()+1;
      	      }
     	      MPI_Bcast(&stringSize, 1, MPI_INT, root, MPI_COMM_WORLD);

      	      char sfilename_c[stringSize];
      	      if(mpiRank == root){
          	  strcpy(sfilename_c,velocityStatsFilename.c_str());
      	      }
 
      	      MPI_Bcast(sfilename_c, stringSize, MPI_CHAR, root, MPI_COMM_WORLD);
      	      velocityStatsFilename.assign(sfilename_c, stringSize);
	  }
      }

      //Same stuff but for thermodynamic stats
      MPI_Bcast(&thermoStats, 1, MPI_C_BOOL, root, MPI_COMM_WORLD);
      if(thermoStats){

	  //if we're doing thermo stats, see if we need to pull from restartt
	  MPI_Bcast(&thermoStatsFromRestart, 1, MPI_C_BOOL, root, MPI_COMM_WORLD);

	  //If we're pulling from restart, get the filename
	  if(thermoStatsFromRestart){

	      if(mpiRank == root){
          	  stringSize = thermoStatsFilename.size()+1;
      	      }
     	      MPI_Bcast(&stringSize, 1, MPI_INT, root, MPI_COMM_WORLD);

      	      char sfilename_c[stringSize];
      	      if(mpiRank == root){
          	  strcpy(sfilename_c, thermoStatsFilename.c_str());
      	      }
 
      	      MPI_Bcast(sfilename_c, stringSize, MPI_CHAR, root, MPI_COMM_WORLD);
      	      thermoStatsFilename.assign(sfilename_c, stringSize);
	  }
      }

      //Type of statistical avering type
      MPI_Bcast(&statsAvgType, 1, MPI_INT, root, MPI_COMM_WORLD);

      //Stats interval
      MPI_Bcast(&stats_interval, 1, MPI_INT, root, MPI_COMM_WORLD);

      //LES Model type
      MPI_Bcast(&lesModel, 1, MPI_INT, root, MPI_COMM_WORLD);
      //LES Model Averaging Type
      MPI_Bcast(&musgsAvgType, 1, MPI_INT, root, MPI_COMM_WORLD);

      //LES DSM Options
      MPI_Bcast(&testFilterType, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&useTaukk, 1, MPI_C_BOOL, root, MPI_COMM_WORLD);
      MPI_Bcast(&dumpCoeffRange, 1, MPI_C_BOOL, root, MPI_COMM_WORLD); 

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

    void parseLESModelString(string vmKey, string inString, LESModel &currentType);
    void parseStatsAvgTypeString(string vmKey, string inString, StatsAvgType &currentType);
    void parseFDTypeFromString(string vmKey, string inString, FDType &currentType);
    void parseFilterTypeFromString(string vmKey, string inString, FilterType &currentType);
    void parseRKTypeFromString(string vmKey, string inString, RKType &currentType);
    void parseTSTypeFromString(string vmKey, string inString, TimeSteppingType &currentType);
    void parseBCTypeFromString(string vmKey, string inString, BCType &currentType);
    void parseBCKindFromString(string vmKey, string inString, BCKind &currentType);
    void bcValidation();

    void parseSpongeFromString(string vmKey, string inString, SpongeKind &spongeKind);

};

#endif 
