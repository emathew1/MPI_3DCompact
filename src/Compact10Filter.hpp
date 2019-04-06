#ifndef _COMPACT10FILTERH_
#define _COMPACT10FILTERH_

#include <math.h>
#include <cstring>
#include <iostream>
#include "AbstractFilter.hpp"

class Compact10Filter: public AbstractFilter{

    public:

 	double alphaF;	
	double a0_10, a1_10, a2_10, a3_10, a4_10, a5_10;

	//Proper boundary filtering schemes...
	double FB1_a, FB1_b, FB1_c, FB1_d, FB1_e, FB1_f, FB1_g, FB1_h, FB1_i, FB1_j, FB1_k;
	double FB2_a, FB2_b, FB2_c, FB2_d, FB2_e, FB2_f, FB2_g, FB2_h, FB2_i, FB2_j, FB2_k;
	double FB3_a, FB3_b, FB3_c, FB3_d, FB3_e, FB3_f, FB3_g, FB3_h, FB3_i, FB3_j, FB3_k;
	double FB4_a, FB4_b, FB4_c, FB4_d, FB4_e, FB4_f, FB4_g, FB4_h, FB4_i, FB4_j, FB4_k;
	double FB5_a, FB5_b, FB5_c, FB5_d, FB5_e, FB5_f, FB5_g, FB5_h, FB5_i, FB5_j, FB5_k;

	//explicit filtering coefficients at boundaries
	double a00, a01, a02, a03, a04;
	double b00, b01, b02, b03, b04;
	double c00, c01, c02, c03, c04;


    Compact10Filter(double alphaF, Domain *dom, BC *bc, Options::BCType bcType, AbstractDerivatives::Direct currentDir){

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

	a0_10    = (193.0 + 126.0*alphaF)/256.0;
	a1_10    = (105.0 + 302.0*alphaF)/256.0;
	a2_10    = (-15.0 +  30.0*alphaF)/64.0;
	a3_10    = ( 45.0 -  90.0*alphaF)/512.0;
	a4_10    = ( -5.0 +  10.0*alphaF)/256.0;
	a5_10    = (  1.0 -   2.0*alphaF)/512.0;
/*
	//F10_1
	FB1_a    = (1023.0 +       alphaF)/1024.0;
	FB1_b    = (   5.0 + 507.0*alphaF)/512.0;
	FB1_c    = ( -45.0 +  45.0*alphaF)/1024.0;
	FB1_d    = (  15.0 -  15.0*alphaF)/128.0;
	FB1_e    = (-105.0 + 105.0*alphaF)/512.0;
	FB1_f    = (  63.0 -  63.0*alphaF)/256.0;
	FB1_g    = (-105.0 + 105.0*alphaF)/512.0;
	FB1_h    = (  15.0 -  15.0*alphaF)/128.0;
	FB1_i    = ( -45.0 +  45.0*alphaF)/1024.0;
	FB1_j    = (   5.0 -   5.0*alphaF)/512.0;
	FB1_k    = (  -1.0 +       alphaF)/1024.0;
*/

	//F10_2
	FB2_a    = (   1.0 +1022.0*alphaF)/1024.0;
	FB2_b    = ( 507.0 +  10.0*alphaF)/512.0;
	FB2_c    = (  45.0 + 934.0*alphaF)/1024.0;
	FB2_d    = ( -15.0 +  30.0*alphaF)/128.0;
	FB2_e    = ( 105.0 - 210.0*alphaF)/512.0;
	FB2_f    = ( -63.0 + 126.0*alphaF)/256.0;
	FB2_g    = ( 105.0 - 210.0*alphaF)/512.0;
	FB2_h    = ( -15.0 +  30.0*alphaF)/128.0;
	FB2_i    = (  45.0 -  90.0*alphaF)/1024.0;
	FB2_j    = (  -5.0 +  10.0*alphaF)/512.0;
	FB2_k    = (   1.0 -   2.0*alphaF)/1024.0;

/*
	//F6_2
	FB2_a    = ( 1.0/64.0  + 31.0*alphaF/32.0);
        FB2_b    = (29.0/32.0  +  3.0*alphaF/16.0);	
	FB2_c    = (15.0/64.0  + 17.0*alphaF/32.0);
	FB2_d    = (-5.0/16.0  +  5.0*alphaF/8.0);
	FB2_e	 = (15.0/64.0  - 15.0*alphaF/32.0);
	FB2_f	 = (-3.0/32.0  +  3.0*alphaF/16.0);
	FB2_g	 = ( 1.0/64.0  -  1.0*alphaF/32.0);
	FB2_h    = 0.0;
	FB2_i    = 0.0;
	FB2_j    = 0.0;
	FB2_k    = 0.0;
*/
/*
	//Interior 2
	FB2_a	 = 0.5*(0.5 + alphaF);
	FB2_b	 =     (0.5 + alphaF);
	FB2_c	 = 0.5*(0.5 + alphaF);
	FB2_d	 = 0.0;
	FB2_e	 = 0.0;
	FB2_f	 = 0.0;
	FB2_g	 = 0.0;
	FB2_h	 = 0.0;
	FB2_i	 = 0.0;
	FB2_j	 = 0.0;
	FB2_k	 = 0.0;
*/

	//F10_3
	FB3_a    = (  -1.0 +   2.0*alphaF)/1024.0;
	FB3_b    = (   5.0 + 502.0*alphaF)/512.0;
	FB3_c    = ( 979.0 +  90.0*alphaF)/1024.0;
	FB3_d    = (  15.0 +  98.0*alphaF)/128.0;
	FB3_e    = (-105.0 + 210.0*alphaF)/512.0;
	FB3_f    = (  63.0 - 126.0*alphaF)/256.0;
	FB3_g    = (-105.0 + 210.0*alphaF)/512.0;
	FB3_h    = (  15.0 -  30.0*alphaF)/128.0;
	FB3_i    = ( -45.0 +  90.0*alphaF)/1024.0;
	FB3_j    = (   5.0 -  10.0*alphaF)/512.0;
	FB3_k    = (  -1.0 +   2.0*alphaF)/1024.0;


/*
	//F8_3
	FB3_a    = (-1.0/256.0 +   1.0*alphaF/128.0);
	FB3_b    = ( 1.0/32.0  +  15.0*alphaF/16.0);
	FB3_c    = (57.0/64.0  +   7.0*alphaF/32.0);
	FB3_d    = ( 7.0/32.0  +   9.0*alphaF/16.0);
	FB3_e 	 =(-35.0/128.0 +  10.0*alphaF/128.0);
	FB3_f	 = ( 7.0/32.0  -   7.0*alphaF/16.0);
	FB3_g	 = (-7.0/64.0  +   7.0*alphaF/32.0);
	FB3_h	 = ( 1.0/32.0  -   1.0*alphaF/16.0);
	FB3_i	 = (-1.0/256.0 +   1.0*alphaF/128.0);
	FB3_j	 = 0.0;
	FB3_k	 = 0.0;
*/
/*
	//F6_3
	FB3_a	 = (-1.0/64.0  +  1.0*alphaF/32.0);
	FB3_b	 = ( 3.0/32.0  + 13.0*alphaF/16.0);
	FB3_c	 = (49.0/64.0  + 15.0*alphaF/32.0);
	FB3_d	 = ( 5.0/16.0  +  3.0*alphaF/8.0);
	FB3_e	 =(-15.0/64.0  + 15.0*alphaF/32.0);
	FB3_f	 = ( 3.0/32.0  -  3.0*alphaF/16.0);
	FB3_g	 = (-1.0/64.0  +  1.0*alphaF/32.0);	
	FB3_h	 = 0.0;
	FB3_i	 = 0.0;
	FB3_j	 = 0.0;
	FB3_k	 = 0.0;
*/
/*
	//Interior 4
	FB3_a	= 0.5*(-1.0/8.0 + 1.0*alphaF/4.0);
	FB3_b	= 0.5*( 1.0/2.0 + 1.0*alphaF/1.0);
	FB3_c	= 5.0/8.0 + 3.0*alphaF/4.0;
	FB3_d	= 0.5*( 1.0/2.0 + 1.0*alphaF/1.0);
	FB3_e	= 0.5*(-1.0/8.0 + 1.0*alphaF/4.0);
	FB3_f	= 0.0;
	FB3_g	= 0.0;
	FB3_h	= 0.0;
	FB3_i	= 0.0;
	FB3_j	= 0.0;
	FB3_k	= 0.0;
*/

	//F10_4
	FB4_a    = (   1.0 -   2.0*alphaF)/1024.0;
	FB4_b    = (  -5.0 +  10.0*alphaF)/512.0;
	FB4_c    = (  45.0 + 934.0*alphaF)/1024.0;
	FB4_d    = ( 113.0 +  30.0*alphaF)/128.0;
	FB4_e    = ( 105.0 + 302.0*alphaF)/512.0;
	FB4_f    = ( -63.0 + 126.0*alphaF)/256.0;
	FB4_g    = ( 105.0 - 210.0*alphaF)/512.0;
	FB4_h    = ( -15.0 +  30.0*alphaF)/128.0;
	FB4_i    = (  45.0 -  90.0*alphaF)/1024.0;
	FB4_j    = (  -5.0 +  10.0*alphaF)/512.0;
	FB4_k    = (   1.0 -   2.0*alphaF)/1024.0;

/*
	//Interior 6
	FB4_a	 = 0.5*( 1.0/32.0 -  1.0*alphaF/16.0);
	FB4_b	 = 0.5*(-3.0/16.0 +  3.0*alphaF/8.0);
	FB4_c	 = 0.5*(15.0/32.0 + 17.0*alphaF/16.0);
	FB4_d	 = 11.0/16.0 + 5.0*alphaF/8.0;
	FB4_e	 = 0.5*(15.0/32.0 + 17.0*alphaF/16.0);
	FB4_f	 = 0.5*(-3.0/16.0 +  3.0*alphaF/8.0);
	FB4_g	 = 0.5*( 1.0/32.0 -  1.0*alphaF/16.0);
	FB4_h    = 0.0;
	FB4_i    = 0.0;
	FB4_j    = 0.0;
	FB4_k    = 0.0;
*/

	//F10_5
	FB5_a    = (  -1.0 +   2.0*alphaF)/1024.0;
	FB5_b    = (   5.0 -  10.0*alphaF)/512.0;
	FB5_c    = ( -45.0 +  90.0*alphaF)/1024.0;
	FB5_d    = (  15.0 +  98.0*alphaF)/128.0;
	FB5_e    = ( 407.0 + 210.0*alphaF)/512.0;
	FB5_f    = (  63.0 + 130.0*alphaF)/256.0;
	FB5_g    = (-105.0 + 210.0*alphaF)/512.0;
	FB5_h    = (  15.0 -  30.0*alphaF)/128.0;
	FB5_i    = ( -45.0 +  90.0*alphaF)/1024.0;
	FB5_j    = (   5.0 -  10.0*alphaF)/512.0;
	FB5_k    = (  -1.0 +   2.0*alphaF)/1024.0;

/*
	//Interior 8
	FB5_a	 = 0.5*(-1.0/128.0 +  1.0*alphaF/64.0);
	FB5_b	 = 0.5*( 1.0/16.0  -  1.0*alphaF/8.0);
	FB5_c	 = 0.5*(-7.0/32.0  + 14.0*alphaF/32.0);
	FB5_d	 = 0.5*( 7.0/16.0  + 18.0*alphaF/16.0);
	FB5_e	 = 93.0/128.0 + 70.0*alphaF/128.0;
	FB5_f	 = 0.5*( 7.0/16.0  + 18.0*alphaF/16.0);
	FB5_g	 = 0.5*(-7.0/32.0  + 14.0*alphaF/32.0);
	FB5_h	 = 0.5*( 1.0/16.0  -  1.0*alphaF/8.0);
	FB5_i	 = 0.5*(-1.0/128.0 +  1.0*alphaF/64.0);
	FB5_j    = 0.0;
	FB5_k    = 0.0;
*/
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
	    //explicitly filter the boundary points (Lele)
	
	    offlowerF[0] = 0.0;
	    //offlowerF[1] = 0.0;
	    //offlowerF[2] = 0.0;
	    offupperF[0] = 0.0;
	    //offupperF[1] = 0.0;
	    //offupperF[2] = 0.0;

	    offlowerF[N-1] = 0.0;
	    //offlowerF[N-2] = 0.0;
	    //offlowerF[N-3] = 0.0;
	    offupperF[N-1] = 0.0;
	    //offupperF[N-2] = 0.0;
	    //offupperF[N-3] = 0.0;

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
