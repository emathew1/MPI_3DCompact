#ifndef _CABSTRACTSOLVERH_
#define _CABSTRACTSOLVERH_

#include "C2Decomp.hpp"
#include "Derivatives.hpp"
#include "Filter.hpp"
#include "IdealGas.hpp"
#include "TimeStepping.hpp"
#include "AbstractSingleBlockMesh.hpp"

//cyclic dependency here...
class AbstractSingleBlockMesh;

class AbstractCSolver{

    public:

	int mpiRank;
	
	int rkStep;
	bool rkLast;
	bool endFlag;

	int timeStep;

	int baseDirection;

	//Stuff for timing functions
	bool useTiming;
	double ft1, ft2;
	 
	//Initial Conditions  
	double *rho0, *p0, *U0, *V0, *W0;

	//Common Solver Data
        double *rho1,  *rhok,  *rhok2,  *rho2;
        double *rhoU1, *rhoUk, *rhoUk2, *rhoU2;
        double *rhoV1, *rhoVk, *rhoVk2, *rhoV2;
        double *rhoW1, *rhoWk, *rhoWk2, *rhoW2;
        double *rhoE1, *rhoEk, *rhoEk2, *rhoE2;       
 
	C2Decomp *c2d;
	Domain *dom;
	BC *bc;
	TimeStepping *ts;
	IdealGas *ig;
        Derivatives *derivX, *derivY, *derivZ;
        Filter *filtX, *filtY, *filtZ;

	//Common Msh Stuff if needed
	AbstractSingleBlockMesh *msh;
	double *J, *J11, *J12, *J13, *J21, *J22, *J23, *J31, *J32, *J33;

	//Each Solver Class needs to have these functions to overwrite the pure virtual ones
	virtual void initializeSolverData() = 0;
	virtual void setInitialConditions() = 0;
	virtual void preStep() = 0;
	virtual void preSubStep() = 0;
	virtual void solveEqnSet() = 0;
	virtual void postSubStep() = 0;
	virtual void updateData() = 0;
	virtual void postStep() = 0;
	virtual void reportAll() = 0;

	//Hook functions
	virtual void subStepTemporalHook(){};
	virtual void fullStepTemporalHook(){};
	virtual void preStepBoundaryHook(){};
	virtual void postStepBoundaryHook(){};

	//RHS source terms
	virtual double contRHSSource(int ip) = 0;
	virtual double xmomRHSSource(int ip) = 0;
	virtual double ymomRHSSource(int ip) = 0;
	virtual double zmomRHSSource(int ip) = 0;
	virtual double engyRHSSource(int ip) = 0;
};

#endif
