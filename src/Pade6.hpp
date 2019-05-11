#ifndef _PADE6H_
#define _PADE6H_

#include <math.h>
#include <cstring>
#include <iostream>
#include "AbstractDerivatives.hpp"

class Pade6: public AbstractDerivatives{


    public:

        //First Derivative

        // Interior Coefficients	
	double alpha_1D;
	double a_1D, b_1D;

	//T6 Dirichlet Coefficients w/ T4 at edge
	double alpha11_1D, alpha21_1D, alpha22_1D;
	double a1_1D, b1_1D, c1_1D, d1_1D, e1_1D, f1_1D;
	double a2_1D, b2_1D, c2_1D, d2_1D, e2_1D; 

        //Second Derivative

        // Interior Coefficients	
	double alpha_2D;
	double a_2D, b_2D;

	//T6 Dirichlet Coefficients w/ T4 at edge
	double alpha11_2D, alpha21_2D, alpha22_2D;
	double a1_2D, b1_2D, c1_2D, d1_2D, e1_2D;
	double a2_2D, b2_2D, c2_2D, d2_2D, e2_2D, f2_2D; 

	int boundaryScheme;

	//Constructor
        Pade6(Domain *dom, Options::BCType bcType, Direct currentDir){

	    //For this method
	    rhsBandwidth = BW5;
	    lhsBandwidth = BW3;

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

	    //Note: Apppears to be stable with 3-4-6-4-3 boundary scheme & 6-6-6-8-10 Filter Scheme for wall bounded
	    //      flows on more complex geometries

	    //1st Derivative coefficients
	    alpha_1D = 1.0/3.0;
	    a_1D     = 14.0/9.0;
	    b_1D     = 1.0/9.0;


	    boundaryScheme = 1;

	    if(boundaryScheme == 1){//Lower-Order @ Node 1

		//Node 1 formulations

	        //3rd @ Node 1
	        alpha11_1D = 2.0;
	        a1_1D = -5.0/2.0;
	        b1_1D =  2.0;
	        c1_1D =  1.0/2.0;
	        d1_1D =  0.0;
	        e1_1D =  0.0;
	        f1_1D =  0.0;

	        /*
	        //Visbal C4 @ Node 1
	        alpha11_1D = 3.0;
	        a1_1D = -17.0/6.0;
	        b1_1D =  3.0/2.0;
	        c1_1D =  3.0/2.0;
	        d1_1D = -1.0/6.0;
	        e1_1D = 0.0;
	        f1_1D = 0.0;
		*/

	 	//Node 2 Formulations

	        //Visbal BC5 @ Node 2
	        alpha21_1D = 1.0/6.0;
	        alpha22_1D = 1.0/2.0;
	        a2_1D = -5.0/9.0;
	        b2_1D = -1.0/2.0;
	        c2_1D =  1.0;
	        d2_1D =  1.0/18.0;

	        /*
	        //Visbal AC5
	        alpha21_1D = 3.0/14.0;
	        alpha22_1D = 3.0/14.0;
	        a2_1D = -19.0/28.0;
	        b2_1D = -5.0/42.0;
	        c2_1D = 6.0/7.0;
	        d2_1D = -1.0/14.0;
	        e2_1D = 1.0/84.0; 
		*/


	   	/*	    
	        //Visbal CC5
	        alpha21_1D = 0.0;
	        alpha22_1D = 3.0/2.0;
	        a2_1D = -1.0/8.0;
	        b2_1D = -11.0/6.0;
	        c2_1D = 3.0/2.0;
	        d2_1D = 1.0/2.0;
	        e2_1D = -1.0/24.0;
		*/


	    }else if(boundaryScheme == 2){//Full Sixth-Order Boundary Treatment

	        alpha11_1D = 5.0;
	        a1_1D = -197.0/60.0;
	        b1_1D =   -5.0/12.0;
	        c1_1D =    5.0;
	        d1_1D =   -5.0/3.0;
 	        e1_1D =    5.0/12.0;
	        f1_1D =   -1.0/20.0;

	        alpha21_1D = 1.0/8.0;
	        alpha22_1D = 3.0/4.0;
	        a2_1D = -43.0/96.0;
	        b2_1D = -5.0/6.0;
	        c2_1D =  9.0/8.0;
	        d2_1D =  1.0/6.0;
	        e2_1D = -1.0/96.0;


	    }else{
		cout << "UNKNOWN BOUNDARY SCHEME PADE6" << endl;
	    }


/*
*/
/*
*/

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
