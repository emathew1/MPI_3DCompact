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
	    //Interior Pentadiagonal 10th order scheme
	    beta  =   1.0/20.0;
	    alpha =   1.0/2.0;
	    a     =  17.0/12.0;
	    b     = 101.0/150.0;
	    c     =   1.0/100.0;

	    //Boundary point 1, using 6th order tri scheme,
	    //can't find stable/correct tenth order for all boundary points
   	    beta1  = 0.0;
   	    alpha1 = 5.0;
	    a1     = -197.0/60.0;
	    b1     =   -5.0/12.0;
	    c1     =    5.0/1.0;
	    d1     =   -5.0/3.0;
 	    e1     =    5.0/12.0;
	    f1     =   -1.0/20.0;
	    g1     =    0.0;
	    h1     =    0.0; 
	    i1     =    0.0; 

	    //Boundary point 2, using 8th order tridiagonal boundary
	    //point from Visbal et al.
	    beta2    = 0.0;
	    alpha2_1 = 1.0/12.0;
	    alpha2_2 = 5.0/4.0;
	    a2       =  -79.0/240.0;
	    b2       =  -77.0/60.0;
	    c2       =   55.0/48.0;
	    d2       =    5.0/9.0;
	    e2       =   -5.0/48.0;
	    f2       =    1.0/60.0;
	    g2       =   -1.0/720.0;
	    h2       =    0.0;

	    //Boundary point 3, using 8th order pentadiagonal interior scheme
	    beta3_1  = 1.0/36.0;
	    beta3_2  = 1.0/36.0;
	    alpha3_1 = 4.0/9.0;
	    alpha3_2 = 4.0/9.0;
	    a3       =  -25.0/54.0/4.0;
	    b3       =  -40.0/27.0/2.0;
	    c3       =    0.0;
	    d3       =   40.0/27.0/2.0;
	    e3       =   25.0/54.0/4.0;
	    f3       =    0.0;
	    g3       =    0.0;


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

    	    //Upper left of the matrix changes
  	    offupper_1D[0] = alpha1;  
	    offupper_1D[1] = alpha2_2;
	    offupper_1D[2] = alpha3_2;

	    offlower_1D[0] = alpha2_1;
	    offlower_1D[1] = alpha3_1;

	    offupper2_1D[0] = beta1;
	    offupper2_1D[1] = beta2;
	    offupper2_1D[2] = beta3_2;

	    offlower2_1D[0] = beta3_1;

	    //Bottom right of the matrix changes
	    offupper_1D[N-2] = alpha2_1;
	    offupper_1D[N-3] = alpha3_1;

	    offlower_1D[N-2] = alpha1;
	    offlower_1D[N-3] = alpha2_2;
	    offlower_1D[N-4] = alpha3_2;

	    offupper2_1D[N-3] = beta3_1;

	    offlower2_1D[N-3] = beta1;
	    offlower2_1D[N-4] = beta2;
	    offlower2_1D[N-5] = beta3_2;
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
