#ifndef _CABSTRACTSINGLEBLOCKMESHH_
#define _CABSTRACTSINGLEBLOCKMESHH_

#include "Derivatives.hpp"

class AbstractSingleBlockMesh{

    public:
	
	int mpiRank;
	double *x, *y, *z;

	double *J;
	double *J11, *J12, *J13;
	double *J21, *J22, *J23;
	double *J31, *J32, *J33;

	double max_xi, max_eta, max_zta;
	double Nx, Ny, Nz;

	Derivatives *derivX, *derivY, *derivZ;

	//Each RK Class needs to have these functions to overwrite the pure virtual ones
	virtual void solveForJacobians() = 0;

};

#endif
