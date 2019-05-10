#ifndef _COMMN5EXPTESTFILTERH_
#define _COMMN5EXPTESTFILTERH_

#include <math.h>
#include <cstring>
#include <iostream>
#include "AbstractFilter.hpp"

//Minimally constrained discrete filters from "COMMUTATIVE FILTERS FOR LES IN COMPLEX GEOMETRIES" (Vasilyev, 1999)
//Not really sure about the wavenumber cutoff for the boundary schemes here

class CommN5ExpTestFilter: public AbstractFilter{

    public:

	double w0, w1, w2, w3; 

	//Boundary points...
	double w0_n1,  w1_n1,  w2_n1, w3_n1, w4_n1, w5_n1; 
	double wn1_n2, w0_n2,  w1_n2, w2_n2, w3_n2, w4_n2; 
	double wn2_n3, wn1_n3, w0_n3, w1_n3, w2_n3, w3_n3; 


    CommN5ExpTestFilter(Domain *dom, BC *bc, Options::BCType bcType, AbstractDerivatives::Direct currentDir){

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

    }

    void FilterPeriodic(double *phi, double *phiF);
    void FilterDirichlet(double *phi, double *phiF);

    void Filter(double *phi, double *phiF);
    void filterField(double *dataIn, double *dataOut);

};

#endif
