#ifndef _PENTA10H_
#define _PENTA10H_

#include <math.h>
#include <cstring>
#include <iostream>
#include "AbstractDerivatives.hpp"

class Penta10: public AbstractDerivatives{


    public:

        //First Derivative

        // Interior Coefficients	
	double alpha, beta;
	double a, b, c;

	//T10 Dirichlet Coefficients

	//Node 1
	double alpha1, beta1;
	double a1, b1, c1, d1, e1, f1, g1, h1, i1;

	//Node 2
	double alpha2_1, alpha2_2, beta2;
	double a2, b2, c2, d2, e2, f2, g2, h2;

	//Node 3
	double alpha3_1, alpha3_2, beta3_1, beta3_2;
	double a3, b3, c3, d3, e3, f3, g3;

	//cp vector for pentadiagonal cyclic solver
	double cpvec[6];

	//Constructor
        Penta10(Domain *dom, Options::BCType bcType, Direct currentDir){

	    //For this method
	    rhsBandwidth = BW7;
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
	    beta  =   1.0/20.0;
	    alpha =   1.0/2.0;
	    a     =  17.0/12.0;
	    b     = 101.0/150.0;
	    c     =   1.0/100.0;

   	    beta1  = 28.0;
   	    alpha1 = 16.0;
	    a1     = -1181.0/280.0;
	    b1     =  -892.0/35.0;
	    c1     =    77.0/5.0;
	    d1     =    56.0/3.0;
 	    e1     =   -35.0/6.0;
	    f1     =    28.0/15.0;
	    g1     =    -7.0/15.0;
	    h1     =     8.0/105.0; 
	    i1     =    -1.0/168.0; 


	    beta2    = 5.0/3.0;
	    alpha2_1 = 1.0/21.0;
	    alpha2_2 = 3.0;
	    a2       = -544.0/2581.0;
	    b2       =  -39.0/20.0;
	    c2       =  -17.0/20.0;
	    d2       =   95.0/36.0;
	    e2       =    5.0/12.0;
	    f2       =   -1.0/20.0;
	    g2       =    1.0/180.0;
	    h2       =   -1.0/2940.0;

	    beta3_1  = 1.0/90.0;
	    beta3_2  = 1.0;
	    alpha3_1 = 4.0/15.0;
	    alpha3_2 = 8.0/9.0;
	    a3       =  -34.0/675.0;
	    b3       = -127.0/225.0;
	    c3       =   -7.0/12.0;
	    d3       =   20.0/27.0;
	    e3       =    4.0/9.0;
	    f3       =    1.0/75.0;
	    g3       =   -1.0/2700.0;


	    diag_1D      = new double[N]; 
	    offlower_1D  = new double[N];
	    offlower2_1D = new double[N];
	    offupper_1D  = new double[N];
	    offupper2_1D = new double[N];


	//Base banded matrix layout
	for(int ip = 0; ip < N; ip++){ 
	    diag_1D[ip] = 1.0;
	    offlower_1D[ip]  = alpha;
	    offupper_1D[ip]  = alpha;
	    offlower2_1D[ip] = beta;
	    offupper2_1D[ip] = beta;
	}

	//These are the modified lhs band values for dirichlet BC's
	if(bcType == Options::DIRICHLET_SOLVE){
	    //NOTE: LAYOUT FOR PENTA SOLVER IS DIFFERENT FROM TRI SOLVER,
	    //TRI SOLVER SKIPS FIRST INDEX FOR LOWER DIAGONAL, PENTA SOLVER 
	    //SKIPS LAST INDEX FOR LOWER DIAGONAL
	    offupper_1D[0] = alpha1;  
	    offupper_1D[1] = alpha2_2;
	    offupper_1D[2] = alpha3_2;

	    offlower_1D[0] = alpha2_1;
	    offlower_1D[1] = alpha3_1;

	    offupper2_1D[0] = beta1;
	    offupper2_1D[1] = beta2;
	    offupper2_1D[2] = beta3_2;

	    offlower2_1D[0] = beta3_1;
	}	

	//cyclical Pentadiagonal solver corner values
	cpvec[0] = beta;
	cpvec[1] = alpha;
	cpvec[2] = beta;
	cpvec[3] = beta;
	cpvec[4] = alpha;
	cpvec[5] = beta;

    }

    //Function's to call...
    void calc1stDerivField(double *dataIn, double *dataOut);
    void calc1stDerivField_TPB(double *dataIn, double *dataOut, double *Nm2, double *Nm1, double *Np1, double *Np2){}; //Empty call since rhs BW is 7
    void calc1stDerivField_TPB(double *dataIn, double *dataOut, double *Nm3, double *Nm2, double *Nm1, double *Np1, double *Np2, double *Np3);

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
