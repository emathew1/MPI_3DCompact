#ifndef _CABSTRACTFILTERH_
#define _CABSTRACTFILTERH_

#include <math.h>
#include <cstring>
#include <iostream>
#include "AbstractDerivatives.hpp"
#include "Macros.hpp"
#include "Domain.hpp"
#include "Utils.hpp"
#include "BC.hpp"

class AbstractFilter{

    public:

	int Nx, Ny, Nz, N;

	int pxSize[3], pySize[3], pzSize[3];
	int pxStart[3], pyStart[3], pzStart[3];
	int pxEnd[3], pyEnd[3], pzEnd[3];

	double *diagF, *offlowerF, *offupperF;

	AbstractDerivatives::Direct currentDir;
	BC::BCType bcType;

	//Functions that a filter object have to define
	virtual void filterField(double *dataIn, double *dataOut) = 0;

};

#endif
