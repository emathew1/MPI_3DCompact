#ifndef _CD2H_
#define _CD2H_

#include <math.h>
#include <cstring>
#include <iostream>
#include "AbstractDerivatives.hpp"

class CD2: public AbstractDerivatives{


    public:

        //First Derivative

        // Interior Coefficients - Second Order Central
	double a, b, c;
	
	// Boundary Coefficients - Second Order Upwind
	double aa, bb, cc;

	//Constructor
        CD2(Domain *dom, Options::BCType bcType, Direct currentDir){

	    //For this method, its really neither but giong to just use it
	    rhsBandwidth = BW5;
	    lhsBandwidth = BW5;

   	    this->Nx = dom->gNx;
	    this->Ny = dom->gNy;
	    this->Nz = dom->gNz;
	    this->dx = dom->dx; 
	    this->dy = dom->dy;
	    this->dz = dom->dz;

	    dom->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);

	    this->currentDir = currentDir;

	    this->bcType = bcType;

	    if(currentDir == DIRX){
		N  = Nx;
		dd = dx;
	    }else if(currentDir == DIRY){
		N  = Ny;
		dd = dy;
	    }else if(currentDir == DIRZ){
		N  = Nz;
		dd = dz; 
	    }

	    //1st Derivative coefficients
	    a     =   -1.0/2.0;
	    b     =    0.0;
	    c     =    1.0/2.0;

	    //Boundary point 1,
	    aa     =   -3.0/2.0;
	    bb     =    4.0/2.0;
	    cc     =   -1.0/2.0; 

    }

    //Function's to call...
    void calc1stDerivField(double *dataIn, double *dataOut);
    void calc1stDerivField_TPB(double *dataIn, double *dataOut, double *Nm2, double *Nm1, double *Np1, double *Np2);
    void calc1stDerivField_TPB(double *dataIn, double *dataOut, double *Nm3, double *Nm2, double *Nm1, double *Np1, double *Np2, double *Np3){};

    void calc1stDeriv(double *phi, double *dphi);
    void calc1stDeriv_TPB(double *phi, double *dphi);

    double calcNeumann(double *f);

    //Need a cleaner way of passing these things...
    void multRHS1stDerivPeriodic(double dh, double *phi, int N, double *RHSvec);
    void multRHS1stDerivPeriodic_TPB(double dh, double *phi, int N, double *RHSvec);
    void multRHS1stDerivDirichlet(double dh, double *phi, int N, double *RHSvec);

    void Compact1stPeriodic(double *phi, double *dphidx);
    void Compact1stPeriodic_TPB(double *phi, double *dphidx);
    void Compact1stDirichlet(double *phi, double *dphidx);
};

#endif
