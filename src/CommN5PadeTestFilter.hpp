#ifndef _COMMN5PADETESTFILTERH_
#define _COMMN5PADETESTFILTERH_

#include <math.h>
#include <cstring>
#include <iostream>
#include "AbstractFilter.hpp"

class CommN5PadeTestFilter: public AbstractFilter{

    public:

 	double alphaF;	

	//Boundary points...
	double w0_n1,  w1_n1,  w2_n1, w3_n1, w4_n1, w5_n1; 
	double wn1_n2, w0_n2,  w1_n2, w2_n2, w3_n2, w4_n2; 
	double wn2_n3, wn1_n3, w0_n3, w1_n3, w2_n3, w3_n3; 

	//Interior Points
	double v0, v1, v2;
	double w0, w1, w2, w3;
	double wi0, wi1, wi2, wi3, wi4, wi5;

	double cpvec[6];

    CommN5PadeTestFilter(Domain *dom, BC *bc, Options::BCType bcType, AbstractDerivatives::Direct currentDir){

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
	offlowerF  = new double[N];
	offlowerF2 = new double[N];
	offupperF  = new double[N];
	offupperF2 = new double[N];

	w0_n1 = 31.0/32.0;
	w1_n1 =  5.0/32.0; 	
	w2_n1 = -5.0/16.0;
	w3_n1 =  5.0/16.0;
	w4_n1 = -5.0/32.0;
	w5_n1 =  1.0/32.0;

	wn1_n2 =  1.0/32.0;
	w0_n2  = 27.0/32.0;
	w1_n2  =  5.0/16.0;
	w2_n2  = -5.0/16.0;
	w3_n2  =  5.0/32.0;
	w4_n2  = -1.0/32.0;

	wn2_n3 = -1.0/32.0;
	wn1_n3 =  5.0/32.0;
	w0_n3  = 11.0/16.0;
	w1_n3  =  5.0/16.0;
	w2_n3  = -5.0/32.0;
	w3_n3  =  1.0/32.0;

	w0 =  11.0/16.0;
	w1 =  15.0/64.0;
	w2 =  -3.0/32.0;
	w3 =   1.0/64.0;

	v0 = 7.0/12.0;
	v1 = 0.0/0.0;
	v2 = 5.0/24.0;

	wi0 =   7.0/24.0;
	wi1 = 175.0/768.0;
	wi2 =   5.0/48.0;
	wi3 =  35.0/1536.0;
	wi4 =   0.0/0.0;
	wi5 =  -1.0/1536.0;
 
        for(int ip = 0; ip < N; ip++){
            diagF[ip] = v0;
            offlowerF[ip]  = v1;
	    offlowerF2[ip] = v2;
            offupperF[ip]  = v1;
	    offupperF2[ip] = v2;
        }
        
	if(bcType == Options::DIRICHLET_SOLVE){
	
	    //Node1 - LHS, explicitly filtered
	    //|1 0 0 ...
	    diagF[0]	  = 1.0;
	    offupperF[0]  = 0.0;
	    offupperF2[0] = 0.0;

	    //Node2 - LHS
	    //|0 1 0 0 ...
	    offlowerF[0]  = 0.0;
	    diagF[1]	  = 1.0;
	    offupperF[1]  = 0.0;
	    offupperF2[1] = 0.0;

	    //Node3 - LHS
	    //|0 0 1 0 0 ...
	    offlowerF2[0] = 0.0;
	    offlowerF[1]  = 0.0;
	    diagF[2]	  = 1.0;
	    offupperF[2]  = 0.0;
	    offupperF[2]  = 0.0;
	
	    //Node4 - LHS
	    //|. 0 0 1 0 0 ...
	    offlowerF2[1] = 0.0;
	    offlowerF[2]  = 0.0;
	    diagF[3]	  = 1.0;
	    offupperF[3]  = 0.0;
	    offupperF2[3] = 0.0;

	    //Node5 - LHS
	    //|. . 0 0 1 0 0 ...
	    offlowerF2[2] = 0.0;
	    offlowerF[3]  = 0.0;
	    diagF[4]	  = 1.0;
	    offupperF[4]  = 0.0;
	    offupperF2[4] = 0.0;
	    
	    //Then by node five we can stop explicitly filtering...
	    //|... v2 v1 v0 v1 v2 ...||f'| = |w5 w4 w3 w2 w1 w0 w1 w2 w3 w4 w5 ...||f|

	    //Now doing the otherside
	    //Node N
	    //|...0 0 1|
	    offlowerF2[N-3] = 0.0;
	    offlowerF[N-2]  = 0.0;
	    diagF[N-1]	    = 1.0;

	    //Node N-1
	    //|...0 0 1 0|
	    offlowerF2[N-4] = 0.0;
	    offlowerF[N-3]  = 0.0;
	    diagF[N-2]	    = 1.0;
	    offupperF[N-2]   = 0.0;
	
	    //Node N-2
	    //|...0 0 1 0 0|
	    offlowerF2[N-5] = 0.0;
	    offlowerF[N-4]  = 0.0;
	    diagF[N-3]	    = 1.0;
	    offupperF[N-3]  = 0.0;
	    offupperF2[N-3] = 0.0;

	    //Node N-3
	    //|...0 0 1 0 0 .|
	    offlowerF2[N-6] = 0.0;
	    offlowerF[N-5]  = 0.0;
	    diagF[N-4]	    = 1.0;
	    offupperF[N-4]  = 0.0;
	    offupperF2[N-4] = 0.0;

	    //Node N-4
	    //|...0 0 1 0 0 . .|
	    offlowerF2[N-7] = 0.0;
	    offlowerF[N-6]  = 0.0;
	    diagF[N-5]	    = 1.0;
	    offupperF[N-5]  = 0.0;
	    offupperF2[N-5] = 0.0;
  	}

	cpvec[0] = v2;
	cpvec[1] = v1;
	cpvec[2] = v2;
	cpvec[3] = v2;
	cpvec[4] = v1;
	cpvec[5] = v2;

    }

    void multRHSPeriodicFilter(double *phi, double *RHSvec);
    void multRHSDirichletFilter(double *phi, double *RHSvec);

    void FilterPeriodic(double *phi, double *phiF);
    void FilterDirichlet(double *phi, double *phiF);

    void compactFilter(double *phi, double *phiF);
    void filterField(double *dataIn, double *dataOut);

};

#endif
