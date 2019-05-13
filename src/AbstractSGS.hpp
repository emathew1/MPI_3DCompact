#ifndef _CABSTRACTSGSH_
#define _CABSTRACTSGSH_

#include "Macros.hpp"
#include "Options.hpp"
#include "C2Decomp.hpp"
#include "AbstractCSolver.hpp"
#include "AbstractSingleBlockMesh.hpp"

class AbstractSGS{

    public:

	int mpiRank;

        //Make local copies for macros...
	int Nx, Ny, Nz, N;

        int pxSize[3], pySize[3], pzSize[3];
        int pxStart[3], pyStart[3], pzStart[3];
        int pxEnd[3], pyEnd[3], pzEnd[3];

 	AbstractFilter *filtX, *filtY, *filtZ;  

	AbstractCSolver *cs;	
	Options::StatsAvgType musgsAvgType;

	double *mu_sgs;	

	virtual void getSGSViscosity(double *gradU[3][3], double *rho, double *rhoU, double *rhoV, double *rhoW, double *rhoE) = 0;

	void doAveraging(double *phi, double *phiF, AbstractFilter *fx, AbstractFilter *fy, AbstractFilter *fz);
	void doAveraging(double *phi, double *phiF);
	void doAveraging();
};

#endif 
