#ifndef _CABSTRACTDERIVATIVESH_
#define _CABSTRACTDERIVATIVESH_

#include <math.h>
#include <cstring>
#include <iostream>
#include "Macros.hpp"
#include "Domain.hpp"
#include "Utils.hpp"
#include "BC.hpp"

class AbstractDerivatives{

    public:

	int Nx, Ny, Nz, N;
	int i, j, k;
	double dx, dy, dz, dd;

	int pxSize[3], pySize[3], pzSize[3];
	int pxStart[3], pyStart[3], pzStart[3];
	int pxEnd[3], pyEnd[3], pzEnd[3];

	double *diag_1D, *offlower_1D, *offupper_1D;
	double *diag_2D, *offlower_2D, *offupper_2D;

	enum Direct {DIRX, DIRY, DIRZ};
	Direct currentDir;

	double Nm3val, Nm2val, Nm1val, N0val, Np1val, Np2val, Np3val; 

	//Functions that a derivative object have to define
	virtual void calc1stDerivField(double *dataIn, double *dataOut) = 0;
	virtual void calc1stDerivField_TPB(double *dataIn, double *dataOut) = 0;
	virtual void calc2ndDerivField(double *dataIn, double *dataOut) = 0;
	

};

#endif
