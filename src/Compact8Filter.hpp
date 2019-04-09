#ifndef _COMPACT8FILTERH_
#define _COMPACT8FILTERH_

#include <math.h>
#include <cstring>
#include <iostream>
#include "AbstractFilter.hpp"

class Compact8Filter: public AbstractFilter{

    public:

 	double alphaF;	
	double a0_8, a1_8, a2_8, a3_8, a4_8;
	double a0_6, a1_6, a2_6, a3_6;

	//Proper boundary filtering schemes
	double FB2_a, FB2_b, FB2_c, FB2_d, FB2_e, FB2_f, FB2_g, FB2_h, FB2_i;
	double FB3_a, FB3_b, FB3_c, FB3_d, FB3_e, FB3_f, FB3_g, FB3_h, FB3_i;
	double FB4_a, FB4_b, FB4_c, FB4_d, FB4_e, FB4_f, FB4_g, FB4_h, FB4_i;

	//explicit filtering coefficients at boundaries
	double a00, a01, a02, a03, a04;
	double b00, b01, b02, b03, b04;
	double c00, c01, c02, c03, c04;


    Compact8Filter(double alphaF, Domain *dom, BC *bc, Options::BCType bcType, AbstractDerivatives::Direct currentDir){

	this->alphaF = alphaF;
	this->Nx = dom->gNx;
	this->Ny = dom->gNy;
	this->Nz = dom->gNz;

	dom->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);

	this->bc = bc;

	this->currentDir = currentDir;

	this->bcType = bcType;

	if(currentDir == AbstractDerivatives::DIRX){
	    N = Nx;
	}else if(currentDir == AbstractDerivatives::DIRY){
	    N = Ny;
	}else if(currentDir == AbstractDerivatives::DIRZ){
	    N = Nz;
	}

	diagF = new double[N];
	offlowerF = new double[N];
	offupperF = new double[N];

	a0_8     = (93.0 + 70.0*alphaF)/128.0;
        a1_8     = ( 7.0 + 18.0*alphaF)/16.0;
        a2_8     = (-7.0 + 14.0*alphaF)/32.0;
        a3_8     = ( 1.0 -  2.0*alphaF)/16.0;
        a4_8     = (-1.0 +  2.0*alphaF)/128.0;

	a0_6     = (11.0 + 10.0*alphaF)/16.0;
	a1_6     = (15.0 + 34.0*alphaF)/32.0;
	a2_6     = (-3.0 +  6.0*alphaF)/16.0;
	a3_6	 = ( 1.0 -  2.0*alphaF)/32.0;	

        FB2_a =   1.0/256.0 + 127.0*alphaF/128.0; 
        FB2_b =  31.0/32.0  +   1.0*alphaF/16.0; 
        FB2_c =   7.0/64.0  +  25.0*alphaF/32.0; 
        FB2_d =  -7.0/32.0  +   7.0*alphaF/16.0; 
        FB2_e =  7.0*(5.0   -  10.0*alphaF)/128.0; 
	FB2_f =  -7.0/32.0  +   7.0*alphaF/16.0;
	FB2_g =   7.0/64.0  -   7.0*alphaF/32.0;
	FB2_h =  -1.0/32.0  +   1.0*alphaF/16.0;
	FB2_i =   1.0/256.0 -   1.0*alphaF/128.0;

        FB3_a =  -1.0/256.0 +   1.0*alphaF/128.0; 
        FB3_b =   1.0/32.0  +  15.0*alphaF/16.0; 
        FB3_c =  57.0/64.0  +   7.0*alphaF/32.0; 
        FB3_d =   7.0/32.0  +   9.0*alphaF/16.0; 
        FB3_e =  7.0*(-5.0  +  10.0*alphaF)/128.0; 
	FB3_f =   7.0/32.0  -   7.0*alphaF/16.0;
	FB3_g =  -7.0/64.0  +   7.0*alphaF/32.0;
	FB3_h =   1.0/32.0  -   1.0*alphaF/16.0;
	FB3_i =  -1.0/256.0 +   1.0*alphaF/128.0;

        FB4_a =   1.0/256.0 -   1.0*alphaF/128.0; 
        FB4_b =  -1.0/32.0  +   1.0*alphaF/16.0; 
        FB4_c =   7.0/64.0  +  25.0*alphaF/32.0; 
        FB4_d =  25.0/32.0  +   7.0*alphaF/16.0; 
        FB4_e =  35.0/128.0 +  29.0*alphaF/64.0; 
	FB4_f =  -7.0/32.0  +   7.0*alphaF/16.0;
	FB4_g =   7.0/64.0  -   7.0*alphaF/32.0;
	FB4_h =  -1.0/32.0  +   1.0*alphaF/16.0;
	FB4_i =   1.0/256.0 -   1.0*alphaF/128.0;

	a00 = 15.0/16.0;	    
	a01 =  4.0/16.0;
	a02 = -6.0/16.0;
	a03 =  4.0/16.0;
	a04 = -1.0/16.0;

	b00 =  1.0/16.0;
	b01 =  3.0/4.0;
	b02 =  6.0/16.0;
	b03 = -4.0/16.0;
	b04 =  1.0/16.0;

	c00 = -1.0/16.0;
	c01 =  4.0/16.0;
	c02 =  5.0/8.0;
	c03 =  4.0/16.0;
	c04 = -1.0/16.0;

        for(int ip = 0; ip < N; ip++){
            diagF[ip] = 1.0;
            offlowerF[ip]  = alphaF;
            offupperF[ip]  = alphaF;
        }
        
	if(bcType == Options::DIRICHLET_SOLVE){
	    offlowerF[0] = 0.0;
	    offupperF[0] = 0.0;

	    offlowerF[N-1] = 0.0;
	    offupperF[N-1] = 0.0;
	}

    }

    void multRHSPeriodicFilter(double *phi, double *RHSvec);
    void multRHSDirichletFilter(double *phi, double *RHSvec);

    void FilterPeriodic(double *phi, double *phiF);
    void FilterDirichlet(double *phi, double *phiF);

    void compactFilter(double *phi, double *phiF);
    void filterField(double *dataIn, double *dataOut);

};

#endif
