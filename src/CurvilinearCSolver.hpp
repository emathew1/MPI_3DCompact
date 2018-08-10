#ifndef _CURVILINEARCSOLVERH_
#define _CURVILINEARCSOLVERH_

#include <iostream>
#include <fstream>
#include <sstream>
#include "Macros.hpp"
#include "Utils.hpp"
#include "SpongeBC.hpp"
#include "AbstractCSolver.hpp"
#include "PngWriter.hpp"
#include "Adt.hpp"

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
        int timeStep;
        double time;
	int filterTimeStep;

	//Kill solver condition
        bool done;

	//non-conserved data
	double  *U_xp, *U, *U_zp, *Ucurv,
		*V_xp, *V, *V_zp, *Vcurv,
		*W_xp, *W, *W_zp, *Wcurv,
	        *T_xp, *T, *T_zp, 
		*p_xp, *p, *p_zp,
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

	double *tempX1,  *tempX2,  *tempX3,  *tempX4, *tempX5;
	double *tempX6,  *tempX7,  *tempX8,  *tempX9, *tempX10;

	double *tempY1,  *tempY2,  *tempY3,  *tempY4,  *tempY5;
	double *tempY6,  *tempY7,  *tempY8,  *tempY9,  *tempY10;
	double *tempY11, *tempY12, *tempY13, *tempY14, *tempY15;

	double *tempZ1,  *tempZ2,  *tempZ3,  *tempZ4,  *tempZ5;
	double *tempZ6,  *tempZ7,  *tempZ8,  *tempZ9,  *tempZ10;

	bool spongeFlag;
	SpongeBC *spg; 

	//See if we have to transpose for Neumann BC calculations
	bool neumannLocalX, neumannLocalZ;

	//Moving Wall BC Velocities
	double X0WallV, X0WallW, X1WallV, X1WallW;
	double Y0WallU, Y0WallW, Y1WallU, Y1WallW;
	double Z0WallU, Z0WallV, Z1WallU, Z1WallV;

	//For drawing images
	PngWriter *pngXY;
	PngWriter *pngXZ;
	PngWriter *pngYZ;

	//Alias'd derivative objects
	Derivatives *derivXi1, *derivXi2, *derivXi3;
	Filter *filtXi1, *filtXi2, *filtXi3;

	//Our Adt search object
	Adt<double> *adt;

	//Constructor to use for this class...
	CurvilinearCSolver(C2Decomp *c2d, Domain *dom, BC *bc, TimeStepping *ts, double alphaF, double mu_ref, bool useTiming){

	    //Take in input information and initialize data structures...
	    this->c2d = c2d;
	    this->dom = dom;
	    this->bc = bc;
	    this->ts = ts;
	    this->alphaF = alphaF;
	    this->mu_ref = mu_ref;
	    this->useTiming = useTiming;

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

	    derivXi1 = derivX;
	    derivXi2 = derivY;
	    derivXi3 = derivZ;

	    //Initialize the filters we're going to use for each direction
	    filtX  = new Filter(alphaF, dom, bc->bcXType, Derivatives::DIRX);
	    filtY  = new Filter(alphaF, dom, bc->bcYType, Derivatives::DIRY);
	    filtZ  = new Filter(alphaF, dom, bc->bcZType, Derivatives::DIRZ);

	    filtXi1 = filtX;
	    filtXi2 = filtY;
	    filtXi3 = filtZ;

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

	    int minYPenXSize;
	    int minYPenZSize;

	    MPI_Allreduce(&pySize[0], &minYPenXSize, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); 
	    MPI_Allreduce(&pySize[2], &minYPenZSize, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

	    IF_RANK0 cout << " > Mininum number of points in Xi direction in Y-Pencil is " << minYPenXSize << ",";
	    if(minYPenXSize >= 7){
		IF_RANK0 cout << " don't need to transpose for Neumann BC's" << endl;
		neumannLocalX = true;
	    }else{
 		IF_RANK0 cout << "need to transpose for Neumann BC's" << endl;
 		IF_RANK0 cout << "WARNING CURRENTLY NOT IMPLETMENTED FOR YBASE SOLVER!!!" << endl;
		neumannLocalX = false;
	    }

	    IF_RANK0 cout << " > Mininum number of points in Zeta direction in Y-Pencil is " << minYPenZSize << ",";
	    if(minYPenZSize >= 7){
		IF_RANK0 cout << " don't need to transpose for Neumann BC's" << endl;
		neumannLocalZ = true;
	    }else{
 		IF_RANK0 cout << "need to transpose for Neumann BC's" << endl;
 		IF_RANK0 cout << "WARNING CURRENTLY NOT IMPLETMENTED FOR YBASE SOLVER!!!" << endl;
		neumannLocalZ = false;
	    }


	    //Initialize the ADT object for interpolation and locating points in the grid...
	    
	    int Nlocal = pySize[0]*pySize[1]*pySize[2];
	    double (*boundBoxMin)[3] = new double[Nlocal][3];
	    double (*boundBoxMax)[3] = new double[Nlocal][3];

	    //Get our coordinates in neighbors across partitions using halo updates
	    double *x_halo = NULL;
	    double *y_halo = NULL;
	    double *z_halo = NULL;

	    //Halo supported only for non-major arrays right now, kind of a pain...
	    double *x_temp1, *y_temp1, *z_temp1;
	    c2d->allocX(x_temp1); 
	    c2d->allocX(y_temp1); 
            c2d->allocX(z_temp1);	

	    MPI_Barrier(MPI_COMM_WORLD);
	    IF_RANK0 cout << "here1" << endl;
	

	    //TODO TODO TODO
	    //This can't be done yet since x, y, & z havent' even been allocated yet!!
	    //This should be moved to some other location in the code after the initilization 
	    // of the mesh object!!!
	    // PERHAPS ADT SHOULD BE A PART OF THE MESH OBJECT, DUHH
	    //TODO TODO TODO	    
 
	    //So move to x-pencil then to non y-major y_pencil...
	    c2d->transposeY2X_MajorIndex(msh->x, x_temp1);
	    c2d->transposeY2X_MajorIndex(msh->y, y_temp1);
	    c2d->transposeY2X_MajorIndex(msh->z, z_temp1);

	    MPI_Barrier(MPI_COMM_WORLD);
	    IF_RANK0 cout << "here2" << endl;

	    //Then back over to y-pencil in x-major array...
	    double *x_temp2, *y_temp2, *z_temp2;
	    c2d->allocY(x_temp2); c2d->allocY(y_temp2); c2d->allocY(z_temp2);	

	    c2d->transposeX2Y(x_temp1, x_temp2);
	    c2d->transposeX2Y(y_temp1, y_temp2);
	    c2d->transposeX2Y(z_temp1, z_temp2);

	    delete[] x_temp1;
	    delete[] y_temp1;
	    delete[] z_temp1;

	    c2d->updateHalo(x_temp2, x_halo, 1, 1);
	    c2d->updateHalo(y_temp2, y_halo, 1, 1);
	    c2d->updateHalo(z_temp2, z_halo, 1, 1);

	    delete[] x_temp2;
	    delete[] y_temp2;
	    delete[] z_temp2;

	    //Cycle through the halo arrays of coordinates
	    for(int kp = 0; kp < pySize[2]; kp++){
		for(int jp = 0; jp < pySize[1]; jp++){
		    for(int ip = 0; ip < pxSize[0]; ip++){

			//This is the non-halo array index
			int ii = kp*pySize[1]*pySize[0] + jp*pySize[0] + ip;

			//Cycle through the points making up a "box" of points

			//This is the halo array index for the same point
		        int iih_0_0_0 = (kp+1)*pySize[1]*(pySize[0]+2) + jp*(pySize[0]+2) + ip + 1;
			
			//Halo array index for i, j, k+1
			int iih_0_0_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (jp)*(pySize[0]+2) + ip + 1;

			//Halo array index for i, j+1, k
			int iih_0_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 1;

			//Halo array index for i, j+1, k+1
			int iih_0_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 1;
			
			//Halo array index for i+1, j, k
			int iih_1_0_0 = (kp+1)*pySize[1]*(pySize[0]+2) + jp*(pySize[0]+2) + ip + 2;

			//Halo array index for i+1, j, k+1
			int iih_1_0_1 = (kp+2)*pySize[1]*(pySize[0]+2) + jp*(pySize[0]+2) + ip + 2;

			//Halo array index for i+1, j+1, k
			int iih_1_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 2;
	
	 		//Halo array index for i+1, j+1, k+1
			int iih_1_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 2;

			if(bc->bcXType == BC::PERIODIC_SOLVE){

			}
	

			if(bc->bcYType == BC::PERIODIC_SOLVE){

			}

			if(bc->bcZType == BC::PERIODIC_SOLVE){

			}

		   }
		}
	    }


	    //adt = new Adt<double>(Nboxes, bbmin, bbmax);	

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
//        	dumpSolution();
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

    	    //updateSponge();
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


