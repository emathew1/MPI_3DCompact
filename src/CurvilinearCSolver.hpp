#ifndef _CURVILINEARCSOLVERH_
#define _CURVILINEARCSOLVERH_

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "Macros.hpp"
#include "Options.hpp"
#include "Utils.hpp"
#include "SpongeBC.hpp"
#include "AbstractCSolver.hpp"
#include "PngWriter.hpp"
#include "CurvilinearInterpolator.hpp"
#include "Pade6.hpp"
#include "Penta10.hpp"
#include "CD2.hpp"
#include "Compact8Filter.hpp"
#include "Compact10Filter.hpp"
#include "Stats.hpp"
#include "AbstractSGS.hpp"
#include "VremanSGS.hpp"

class CurvilinearCSolver: public AbstractCSolver{

    public:

	double alphaF;
	double mu_ref;

	//Make local copies for macros...
	int Nx, Ny, Nz, N;

        int pxSize[3], pySize[3], pzSize[3];
        int pxStart[3], pyStart[3], pzStart[3];
        int pxEnd[3], pyEnd[3], pzEnd[3];

	int totRank;

	double t1, t2;

	
	
	//Track the current time and timestep
	int filterTimeStep;

	//Kill solver condition
        bool done;

	//non-conserved data
	double  *U_xp, *U_zp, *Ucurv,
		*V_xp, *V_zp, *Vcurv,
		*W_xp, *W_zp, *Wcurv,
	        *T_xp, *T_zp, 
		*p_xp, *p_zp,
	       *mu,
	      *sos;

	//derivatives of data
	double *dU1, *dU2, *dU3;
	double *dV1, *dV2, *dV3;
	double *dW1, *dW2, *dW3;
	double *dT1, *dT2, *dT3;


	double *Tau11,  *Tau12,  *Tau13;
	double *Tau21,  *Tau22,  *Tau23;
	double *Tau31,  *Tau32,  *Tau33;

	double *cont_1, *cont_2, *cont_3;
	double *mom1_1, *mom1_2, *mom1_3;
	double *mom2_1, *mom2_2, *mom2_3;
	double *mom3_1, *mom3_2, *mom3_3;
	double *engy_1, *engy_2, *engy_3;

	vector<double*> tempXVec;

	vector<double*> tempYVec;

	vector<double*> tempZVec;

	//Stuff for different styles of interprocessor computation...
	enum CompStyle {VANILLA, OCC, VANILLA_CHUNKED, OCC_CHUNKED};
        CompStyle compStyle;

	//Vector containers for send and receive buffers for OCC style
	vector<double*> sbufVec;
	vector<double*> rbufVec;
	
	//Field for visualizing the processor decomposition of the field
	double *rankfield;

	bool spongeFlag;
	SpongeBC *spg; 

	//See if we have to transpose for Neumann BC calculations
	bool neumannLocalX, neumannLocalZ;

	//Moving Wall BC Velocities
	double X0WallV, X0WallW, X1WallV, X1WallW;
	double Y0WallU, Y0WallW, Y1WallU, Y1WallW;
	double Z0WallU, Z0WallV, Z1WallU, Z1WallV;

	//For drawing images
	list<PngWriter*> imageList;

	//Statstics object
	Stats *stats;
	bool statsFlag;

	//LES Object
	AbstractSGS *les;
	bool LESFlag;

	//Alias'd derivative objects
	AbstractDerivatives *derivXi1, *derivXi2, *derivXi3;
	AbstractFilter *filtXi1, *filtXi2, *filtXi3;

	//Constructor to use for this class...
	CurvilinearCSolver(C2Decomp *c2d, Domain *dom, BC *bc, TimeStepping *ts, Options *opt){

	    //Take in input information and initialize data structures...
	    this->c2d = c2d;
	    this->dom = dom;
	    this->bc = bc;
	    this->ts = ts;
	    this->opt = opt;
	    this->alphaF = opt->alphaF;
	    this->mu_ref = opt->mu_ref;
	    this->useTiming = opt->useTiming;

	    baseDirection = 1;

	    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  	    MPI_Comm_size(MPI_COMM_WORLD, &totRank);	

	    ig = new IdealGas(dom, mu_ref);

	    //give ourselves the local copies of the domain sizes
	    Nx = dom->gNx;
	    Ny = dom->gNy;
	    Nz = dom->gNz;
	    N  = dom->gN;
            dom->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);

	    //initialize time and timestep
   	    time = opt->time;
	    timeStep = opt->timeStep;
	    filterTimeStep = 0;
	    endFlag = false;
	    done = false;
	    rkLast = false;

	    //Computational style, should be in inputs
	    compStyle = OCC;
	
	    //Allocate our arrays for the solver data
	    initializeSolverData();		    	    

	    //Initialize our derivative calculations for each direction...

	    if(opt->xFDType == Options::CD2){
	        derivX = new CD2(dom, bc->bcXType, AbstractDerivatives::DIRX);
	    }else if(opt->xFDType == Options::PADE6){
	        derivX = new Pade6(dom, bc->bcXType, AbstractDerivatives::DIRX);
	    }else if(opt->xFDType == Options::PENTA10){
	        derivX = new Penta10(dom, bc->bcXType, AbstractDerivatives::DIRX);
	    }else{
		cout << "Should never get here? unknown x-derivative" << endl;
		MPI_Abort(MPI_COMM_WORLD, -10);
	    }

	    if(opt->yFDType == Options::CD2){
	        derivY = new CD2(dom, bc->bcYType, AbstractDerivatives::DIRY);
	    }else if(opt->yFDType == Options::PADE6){
	        derivY = new Pade6(dom, bc->bcYType, AbstractDerivatives::DIRY);
	    }else if(opt->yFDType == Options::PENTA10){
	        derivY = new Penta10(dom, bc->bcYType, AbstractDerivatives::DIRY);
	    }else{
		cout << "Should never get here? unknown y-derivative" << endl;
		MPI_Abort(MPI_COMM_WORLD, -10);
	    }

	    if(opt->zFDType == Options::CD2){
	        derivZ = new CD2(dom, bc->bcZType, AbstractDerivatives::DIRZ);
	    }else if(opt->zFDType == Options::PADE6){
	        derivZ = new Pade6(dom, bc->bcZType, AbstractDerivatives::DIRZ);
	    }else if(opt->zFDType == Options::PENTA10){
	        derivZ = new Penta10(dom, bc->bcZType, AbstractDerivatives::DIRZ);
	    }else{
		cout << "Should never get here? unknown z-derivative" << endl;
		MPI_Abort(MPI_COMM_WORLD, -10);
	    }

	    derivXi1 = derivX;
	    derivXi2 = derivY;
	    derivXi3 = derivZ;

	    //Initialize the filters we're going to use for each direction
	    if(opt->filterType == Options::COMPACT8){
	        filtX  = new Compact8Filter(alphaF, dom, bc, bc->bcXType, AbstractDerivatives::DIRX);
	        filtY  = new Compact8Filter(alphaF, dom, bc, bc->bcYType, AbstractDerivatives::DIRY);
	        filtZ  = new Compact8Filter(alphaF, dom, bc, bc->bcZType, AbstractDerivatives::DIRZ);
	    }else if(opt->filterType == Options::COMPACT10){
	        filtX  = new Compact10Filter(alphaF, dom, bc, bc->bcXType, AbstractDerivatives::DIRX);
	        filtY  = new Compact10Filter(alphaF, dom, bc, bc->bcYType, AbstractDerivatives::DIRY);
	        filtZ  = new Compact10Filter(alphaF, dom, bc, bc->bcZType, AbstractDerivatives::DIRZ);

	    }else{
		cout << "Should never get here? unknown z-derivative" << endl;
		MPI_Abort(MPI_COMM_WORLD, -10);

	    }

	    filtXi1 = filtX;
	    filtXi2 = filtY;
	    filtXi3 = filtZ;

 	    X0WallV = 0.0; X0WallW = 0.0; X1WallV = 0.0; X1WallW = 0.0;
	    Y0WallU = 0.0; Y0WallW = 0.0; Y1WallU = 0.0; Y1WallW = 0.0;
	    Z0WallU = 0.0; Z0WallV = 0.0; Z1WallU = 0.0; Z1WallV = 0.0;

	    //Initialize the statistics object
	    statsFlag = opt->velocityStats || opt->thermoStats;
	    if(statsFlag){
	        stats = new Stats(this, opt);
	    }else{
		stats = NULL;
	    }

	    //Initialize the LES stuff if we need to...
	    if(opt->lesModel == Options::NOMODEL){
		LESFlag = false;
		les = NULL;
	    }else if(opt->lesModel == Options::VREMAN){
		LESFlag = true;
		les = new VremanSGS(this);
	    } 

	    t1 = MPI_Wtime();
	}

	//Pre solver utility functions
	void addImageOutput(PngWriter *pw);
	void generateImagePlane(PngWriter *pw);

	//Pre Step Functions...
	void calcDtFromCFL();

	//Pre Sub Step Functions...
	void preStepBCHandling();
	void preStepDerivatives();

	//Solve Eqn Set Functions...
	void solveContinuity();
	void solveXMomentum();
	void solveYMomentum();
	void solveZMomentum();
	void solveEnergy();

	//Post Sub Step Functions...
	void postStepBCHandling();
	void filterConservedData();
	void updateNonConservedData();

	//Post Step Functions
	void updateSponge();
	void writeImages();
	void writePlaneImageForVariable(PngWriter *pw);
	void checkSolution();
	void dumpSolution();
	void checkEnd();
	void reportAll();

	//Hook functions
	void initialHook();
	void fullStepTemporalHook();
	void subStepTemporalHook();
	void preStepBoundaryHook();
	void postStepBoundaryHook();

	double contRHSSource(int ip);
	double xmomRHSSource(int ip);
	double ymomRHSSource(int ip);
	double zmomRHSSource(int ip);
	double engyRHSSource(int ip);

	//Inline, Or general solver functions
	inline double calcSpongeSource(double phi, double phiSpongeAvg, double sigma){
        	return sigma*(phiSpongeAvg - phi);
	};

	//Derivative calculation functions
	void computeGradient(vector<double*> vecIn, vector<double*>vecOut);
	void computeGradDotComponents(vector<double*> vecIn, vector<double*>vecOut);

	bool checkForAndDeleteKillFile(string killFileName);

	/////////////////////////////////////
	//Our Generalized Solver Functions //
	/////////////////////////////////////

	void setInitialConditions();
	void initializeSolverData();

	void preStep(){

   	    if(timeStep == opt->timeStep){
//        	dumpSolution();
        	writeImages();
    	    }
    	    calcDtFromCFL();

	}

	void preSubStep(){

    	    preStepBCHandling();
	    preStepBoundaryHook();

    	    preStepDerivatives();

	}

	void solveEqnSet(){

    	    solveContinuity();
            solveXMomentum();
    	    solveYMomentum();
    	    solveZMomentum();
    	    solveEnergy();
	    
	}

	void postSubStep(){

    	    postStepBCHandling();
	    postStepBoundaryHook();
	    subStepTemporalHook();
	}

	void updateData(){
 
     	    if(rkLast || opt->subStepFiltering){
	       	filterConservedData();
   	    }
 
	    updateNonConservedData();
	}

	void postStep(){

	    fullStepTemporalHook();

    	    //updateSponge();
    	    checkSolution();

	    if(timeStep%opt->stats_interval == 0 && statsFlag){
		stats->updateStatsFields();
	    } 


    	    if(timeStep%ts->dumpStep == 0){
        	dumpSolution();
	    }

       	    writeImages();

    	    checkEnd();

	}



};



#endif


