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

	//Constructor
        Penta10(Domain *dom, Options::BCType bcType, Direct currentDir){

	    //For this method
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
	    alpha =   1.0/2.0;
	    beta  =   1.0/20.0;
	    a     =  17.0/12.0;
	    b     = 101.0/150.0;
	    c     =   1.0/100.0;

   	    alpha1 = 16.0;
   	    beta1  = 28.0;
	    a1     = -1181.0/280.0;
	    b1     =  -892.0/35.0;
	    c1     =    77.0/5.0;
	    d1     =    56.0/3.0;
 	    e1     =   -35.0/6.0;
	    f1     =    28.0/15.0;
	    g1     =    -7.0/15.0;
	    h1     =     8.0/105.0; 
	    i1     =    -1.0/168.0; 

	    //LEFT OFF HERE

	    alpha21_1D = 1.0/8.0;
	    alpha22_1D = 3.0/4.0;
	    a2_1D = -43.0/96.0;
	    b2_1D = -5.0/6.0;
	    c2_1D =  9.0/8.0;
	    d2_1D =  1.0/6.0;
	    e2_1D = -1.0/96.0;

	    //2nd Derivative coefficients
	    alpha_2D = 2.0/11.0;
	    a_2D     = 12.0/11.0;
	    b_2D     = 3.0/11.0;

	    //4th order here
   	    alpha11_2D = 10.0;
	    a1_2D =  145.0/12.0;
	    b1_2D =  -76.0/3.0;
	    c1_2D =   29.0/2.0;
	    d1_2D =   -4.0/3.0;
 	    e1_2D =    1.0/12.0;

	    //6th order here...
	    alpha21_2D = 2.0/11.0;
	    alpha22_2D = -131.0/22.0;
	    a2_2D =  177.0/88.0;
	    b2_2D = -507.0/44.0;
	    c2_2D =  783.0/44.0;
	    d2_2D = -201.0/22.0;
	    e2_2D =  81.0/88.0;
	    f2_2D =  -3.0/44.0;


	    diag_1D     = new double[N]; 
	    offlower_1D = new double[N];
	    offupper_1D = new double[N];

	    diag_2D     = new double[N]; 
	    offlower_2D = new double[N];
	    offupper_2D = new double[N];

	//Left off here!     
  
	for(int ip = 0; ip < N; ip++){ 
	    diag_1D[ip] = 1.0;
	    offlower_1D[ip]  = alpha_1D;
	    offupper_1D[ip]  = alpha_1D;
	}

	for(int ip = 0; ip < N; ip++){
	    diag_2D[ip] = 1.0;
	    offlower_2D[ip]  = alpha_2D;
	    offupper_2D[ip]  = alpha_2D;
	}

	if(bcType == Options::DIRICHLET_SOLVE){
	    offupper_1D[0] = alpha11_1D;  
	    offupper_1D[1] = alpha22_1D;
	    offlower_1D[1] = alpha21_1D;

	    offlower_1D[N-1] = alpha11_1D;
	    offlower_1D[N-2] = alpha22_1D;
	    offupper_1D[N-2] = alpha21_1D;

	    offupper_2D[0] = alpha11_2D;  
	    offupper_2D[1] = alpha22_2D;
	    offlower_2D[1] = alpha21_2D;

	    offlower_2D[N-1] = alpha11_2D;
	    offlower_2D[N-2] = alpha22_2D;
	    offupper_2D[N-2] = alpha21_2D;
	}	


    }

    //Function's to call...
    void calc1stDerivField(double *dataIn, double *dataOut);
    void calc1stDerivField_TPB(double *dataIn, double *dataOut, double *Nm2, double *Nm1, double *Np1, double *Np2);
    void calc1stDerivField_TPB(double *dataIn, double *dataOut, double *Nm3, double *Nm2, double *Nm1, double *Np1, double *Np2, double *Np3){};//empty call here since our bandwidth is 5
    void calc2ndDerivField(double *dataIn, double *dataOut);

    void calc1stDeriv(double *phi, double *dphi);
    void calc1stDeriv_TPB(double *phi, double *dphi);
    void calc2ndDeriv(double *phi, double *dphi);

    double calcNeumann(double *f);

    //Need a cleaner way of passing these things...
    void multRHS1stDerivPeriodic(double dh, double *phi, int N, double *RHSvec);
    void multRHS1stDerivPeriodic_TPB(double dh, double *phi, int N, double *RHSvec);
    void multRHS2ndDerivPeriodic(double dh, double *phi, int N, double *RHSvec);
    void multRHS1stDerivDirichlet(double dh, double *phi, int N, double *RHSvec);
    void multRHS2ndDerivDirichlet(double dh, double *phi, int N, double *RHSvec);


    void Compact1stPeriodic(double *phi, double *dphidx);
    void Compact1stPeriodic_TPB(double *phi, double *dphidx);
    void Compact2ndPeriodic(double *phi, double *dphidx);
    void Compact1stDirichlet(double *phi, double *dphidx);
    void Compact2ndDirichlet(double *phi, double *dphidx);
};

#endif
