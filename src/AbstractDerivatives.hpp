#ifndef _CABSTRACTDERIVATIVESH_
#define _CABSTRACTDERIVATIVESH_

#include <math.h>
#include <cstring>
#include <iostream>
#include "Macros.hpp"
#include "Options.hpp"
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
	double *offlower2_1D, *offupper2_1D;
	double *diag_2D, *offlower_2D, *offupper_2D;

	enum Direct {DIRX, DIRY, DIRZ};
	Direct currentDir;

	Options::BCType bcType;

	double Nm3val, Nm2val, Nm1val, N0val, Np1val, Np2val, Np3val; 

	enum BW {BW3, BW5, BW7};
	BW lhsBandwidth;
	BW rhsBandwidth;

	//Functions that a derivative object have to define
	virtual void calc1stDerivField(double *dataIn, double *dataOut) = 0;

	//SLOPPY OVERLOADED IMPLEMENTATION INCOMING
	//transformed periodic boundaries for rhsBandwidth of 5
	virtual void calc1stDerivField_TPB(double *dataIn, double *dataOut, double *Nm2, double *Nm1, double *Np1, double *Np2) = 0;
	//transformed periodic boundaries for rhsBandwidth of 7
	virtual void calc1stDerivField_TPB(double *dataIn, double *dataOut, double *Nm3, double *Nm2, double *Nm1, double *Np1, double *Np2, double *Np3) = 0;

	//Not necessary with the current solver that just does two 1st derivs
	//virtual void calc2ndDerivField(double *dataIn, double *dataOut) = 0;
	
	//Going to move the calc Neumann function here to force it to be defined for each kind of derivative
	virtual double calcNeumann(double *f) = 0;

};

#endif
