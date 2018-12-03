#ifndef _UNIFORMCSOLVERH_
#define _UNIFORMCSOLVERH_

#include <iostream>
#include <fstream>
#include <sstream>
#include "Macros.hpp"
#include "Utils.hpp"
#include "SpongeBC.hpp"
#include "AbstractCSolver.hpp"
#include "PngWriter.hpp"

class UniformCSolver: public AbstractCSolver{

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
        int timeStep;
        double time;
	int filterTimeStep;

	//Kill solver condition
        bool done;

	//non-conserved data
	double  *U, *U2, *U3, 
		*V, *V2, *V3, 
		*W, *W2, *W3,
	        *T, *T2, *T3, 
		*p, *p2, *p3,
	       *mu,
	      *Amu,
		*k, 
	      *sos;

	//derivatives of data
	double *Ux, *Uy, *Uz;
	double *Vx, *Vy, *Vz;
	double *Wx, *Wy, *Wz;

	double *Uxx, *Uyy, *Uzz;
	double *Vxx, *Vyy, *Vzz;
	double *Wxx, *Wyy, *Wzz;

	double *Uxy, *Uxz, *Uyz;
	double *Vxy, *Vxz, *Vyz;
	double *Wxy, *Wxz, *Wyz;

	double *Tx,  *Ty,  *Tz;
	double *Txx, *Tyy, *Tzz;

	double *contEulerX, *contEulerY, *contEulerZ;
	double *momXEulerX, *momXEulerY, *momXEulerZ;
	double *momYEulerX, *momYEulerY, *momYEulerZ;
	double *momZEulerX, *momZEulerY, *momZEulerZ;
	double *engyEulerX, *engyEulerY, *engyEulerZ;

	double *tempX1,  *tempX2,  *tempX3,  *tempX4, *tempX5;

	double *tempY1,  *tempY2,  *tempY3,  *tempY4,  *tempY5;
	double *tempY6,  *tempY7,  *tempY8,  *tempY9,  *tempY10;
	double *tempY11, *tempY12, *tempY13, *tempY14, *tempY15;
	double *tempY16, *tempY17, *tempY18, *tempY19, *tempY20;
	double *tempY21, *tempY22, *tempY23, *tempY24, *tempY25;
	double *tempY26, *tempY27, *tempY28;

	double *tempZ1,  *tempZ2,  *tempZ3,  *tempZ4,  *tempZ5;
	double *tempZ6,  *tempZ7,  *tempZ8,  *tempZ9,  *tempZ10;
	double *tempZ11, *tempZ12, *tempZ13, *tempZ14, *tempZ15;
	double *tempZ16, *tempZ17, *tempZ18, *tempZ19, *tempZ20;
	double *tempZ21, *tempZ22, *tempZ23, *tempZ24, *tempZ25;
	double *tempZ26, *tempZ27, *tempZ28, *tempZ29, *tempZ30;
	double *tempZ31, *tempZ32, *tempZ33, *tempZ34;

	//TESTING
        double *sbuf1, *sbuf2, *sbuf3, *sbuf4, *sbuf5;
	double *rbuf1, *rbuf2, *rbuf3, *rbuf4, *rbuf5;


	double *temp, *temp2, *temp3, *temp4;
	double *transRho, *transRhoU, *transRhoV, *transRhoW, *transRhoE;
	double *transUx, *transVx, *transWx;
	double *transUy, *transVy, *transWy;

	//New memory allocated for AWS solver...
        double *transTempUy;
        double *transTempVy;
    	double *transTempWy;
    	double *transTempUyy;
    	double *transTempVyy;
    	double *transTempWyy;
    	double *transTempTy; 
    	double *transTempTyy; 
	double *transTempUxy;
	double *transTempVxy;
	double *transTempWxy;
	double *transTempUyz;
	double *transTempVyz;
	double *transTempWyz;
	double *transTempContEuler;
	double *transTempXEuler;
	double *transTempYEuler;
	double *transTempZEuler;
	double *transTempEngEuler;

	bool spongeFlag;
	SpongeBC *spg; 

	//See if we have to transpose for Neumann BC calculations
	bool neumannLocalY, neumannLocalZ;

	//Moving Wall BC Velocities
	double X0WallV, X0WallW, X1WallV, X1WallW;
	double Y0WallU, Y0WallW, Y1WallU, Y1WallW;
	double Z0WallU, Z0WallV, Z1WallU, Z1WallV;

	//For drawing images
	PngWriter *pngXY;
	PngWriter *pngXZ;
	PngWriter *pngYZ;

	//Constructor to use for this class...
	UniformCSolver(C2Decomp *c2d, Domain *dom, BC *bc, TimeStepping *ts, double alphaF, double mu_ref, bool useTiming){

	    //Take in input information and initialize data structures...
	    this->c2d = c2d;
	    this->dom = dom;
	    this->bc = bc;
	    this->ts = ts;
	    this->alphaF = alphaF;
	    this->mu_ref = mu_ref;
	    this->useTiming = useTiming;

	    baseDirection = 0;

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
   	    time = 0.0;
	    timeStep = 0;
	    filterTimeStep = 0;
	    endFlag = false;
	    done = false;
	    rkLast = false;

	    //Allocate our arrays for the solver data
	    initializeSolverData();		    	    


	    //Initialize the sponge boundary conditions if necessary
	    if(bc->bcX0 == BC::SPONGE || bc->bcX1 == BC::SPONGE || bc->bcY0 == BC::SPONGE || bc->bcY1 == BC::SPONGE || bc->bcZ0 == BC::SPONGE || bc->bcZ1 == BC::SPONGE){
		spongeFlag = true;
		spg = new SpongeBC(dom, ig, bc, c2d, baseDirection, mpiRank);
	    }else{
		spongeFlag = false;
		spg = NULL;
	    }

	    //Initialize our derivative calculations for each direction...
	    derivX = new Derivatives(dom, bc->bcXType, Derivatives::DIRX);
	    derivY = new Derivatives(dom, bc->bcYType, Derivatives::DIRY);
	    derivZ = new Derivatives(dom, bc->bcZType, Derivatives::DIRZ);

	    //Initialize the filters we're going to use for each direction
	    filtX  = new Filter(alphaF, dom, bc->bcXType, Derivatives::DIRX);
	    filtY  = new Filter(alphaF, dom, bc->bcYType, Derivatives::DIRY);
	    filtZ  = new Filter(alphaF, dom, bc->bcZType, Derivatives::DIRZ);

 	    X0WallV = 0.0; X0WallW = 0.0; X1WallV = 0.0; X1WallW = 0.0;
	    Y0WallU = 0.0; Y0WallW = 0.0; Y1WallU = 0.0; Y1WallW = 0.0;
	    Z0WallU = 0.0; Z0WallV = 0.0; Z1WallU = 0.0; Z1WallV = 0.0;

	    //only rank 0 write images for now...
	    IF_RANK0{
		pngXY = new PngWriter(Nx, Ny); 
		pngXZ = new PngWriter(Nz, Nx); 
		pngYZ = new PngWriter(Ny, Nz); 
	    }else{
		pngXY = NULL;
		pngXZ = NULL;
		pngYZ = NULL;
	    }
	

	    int minXPenYSize;
	    int minXPenZSize;

	    MPI_Allreduce(&pxSize[1], &minXPenYSize, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); 
	    MPI_Allreduce(&pxSize[2], &minXPenZSize, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

	    IF_RANK0 cout << " > Mininum number of points in Y direction in X-Pencil is " << minXPenYSize << ",";
	    if(minXPenYSize >= 7){
		IF_RANK0 cout << " don't need to transpose for Neumann BC's" << endl;
		neumannLocalY = true;
	    }else{
 		IF_RANK0 cout << "need to transpose for Neumann BC's" << endl;
		neumannLocalY = false;
	    }

	    IF_RANK0 cout << " > Mininum number of points in Z direction in X-Pencil is " << minXPenZSize << ",";
	    if(minXPenZSize >= 7){
		IF_RANK0 cout << " don't need to transpose for Neumann BC's" << endl;
		neumannLocalZ = true;
	    }else{
 		IF_RANK0 cout << "need to transpose for Neumann BC's" << endl;
		neumannLocalZ = false;
	    }



	    t1 = MPI_Wtime();
	}

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
	void writePlaneImageForVariable(double *var, string varName, string timeStepString, int plane, double fraction);
	void checkSolution();
	void dumpSolution();
	void checkEnd();
	void reportAll();

	virtual void temporalCalculations(){};

	//Inline, Or general solver functions
	inline double calcSpongeSource(double phi, double phiSpongeAvg, double sigma){
        	return sigma*(phiSpongeAvg - phi);
	};

	/////////////////////////////////////
	//Our Generalized Solver Functions //
	/////////////////////////////////////

	void setInitialConditions();
	void initializeSolverData();

	void preStep(){
   	    if(timeStep == 0){
        	dumpSolution();
        	writeImages();
    	    }
    	    calcDtFromCFL();
	}

	void preSubStep(){
    	    preStepBCHandling();
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
	}

	void updateData(){
    	    if(rkLast){
        	filterConservedData();
    	    }
    	    updateNonConservedData();

	}

	void postStep(){

    	    updateSponge();
    	    checkSolution();

    	    if(timeStep%ts->dumpStep == 0)
        	dumpSolution();

    	    if(timeStep%ts->imageStep == 0)
        	writeImages();

	    temporalCalculations();

    	    checkEnd();
	}



};



#endif


