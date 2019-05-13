#ifndef _CDSMSGSH_
#define _CDSMSGSH_

#include "Macros.hpp"
#include "AbstractFilter.hpp"
#include "AbstractDerivatives.hpp"
#include "CommN5ExpTestFilter.hpp"
#include "CommN5PadeTestFilter.hpp"

class DSMSGS: public AbstractSGS{

  public:

    DSMSGS(AbstractCSolver *cs){
	
	this->mpiRank = cs->mpiRank;
	this->cs      = cs;

	Nx = cs->dom->gNx;
	Ny = cs->dom->gNy;
	Nz = cs->dom->gNz;
	N  = cs->dom->gN;
        cs->dom->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);

	musgsAvgType = cs->opt->musgsAvgType;

	//initialize mu_sgs
	cs->c2d->allocY(mu_sgs);

	//This should be an input option...
	int filterType = 1;

        if(filterType == 1){
	    filtX = new CommN5PadeTestFilter(cs->dom, cs->bc, cs->bc->bcXType, AbstractDerivatives::DIRX);
	    filtY = new CommN5PadeTestFilter(cs->dom, cs->bc, cs->bc->bcYType, AbstractDerivatives::DIRY);
	    filtZ = new CommN5PadeTestFilter(cs->dom, cs->bc, cs->bc->bcZType, AbstractDerivatives::DIRZ);
	}else if(filterType == 2{
	    filtX = new CommN5ExpTestFilter(cs->dom, cs->bc, cs->bc->bcXType, AbstractDerivatives::DIRX);
	    filtY = new CommN5ExpTestFilter(cs->dom, cs->bc, cs->bc->bcYType, AbstractDerivatives::DIRY);
	    filtZ = new CommN5ExpTestFilter(cs->dom, cs->bc, cs->bc->bcZType, AbstractDerivatives::DIRZ);
	}

    }

    void filterQuantity(double *phi, double *phiF){

	fy->filterField(phi, cs->tempY1);
	cs->c2d->transposeY2Z_MajorIndex(cs->tempY1, cs->tempZ1);
	fz->filterField(cs->tempZ1, cs->tempZ2);
	cs->c2d->transposeZ2Y_MajorIndex(cs->tempZ2, cs->tempY1);
	cs->c2d->transposeY2X_MajorIndex(cs->tempY1, cs->tempX1);
	fx->filterField(cs->tempX1, cs->tempX2);
	cs->c2d->transposeX2Y_MajorIndex(cs->tempX2, phiF);
    }
   
    void getSGSViscosity(double *gradU[3][3], double *rho, double *rhoU, double *rhoV, double *rhoW, double *rhoE){


	    double *S00, *S01, *S02;
	    double       *S11, *S12;
	    double 	       *S22; 
	    double *Smag;	    

	    cs->c2d->allocY(S00);
	    cs->c2d->allocY(S01);
	    cs->c2d->allocY(S02);
	    cs->c2d->allocY(S11);
	    cs->c2d->allocY(S12);
	    cs->c2d->allocY(S22);
	    cs->c2d->allocY(Smag);

	    FOR_XYZ_YPEN{
		S00[ip] = gradU[0][0][ip];
		S11[ip] = gradU[1][1][ip];
		S22[ip] = gradU[2][2][ip];
		S01[ip] = 0.5*(gradU[0][1][ip] + gradU[1][0][ip]); 	
		S02[ip] = 0.5*(gradU[0][2][ip] + gradU[2][0][ip]); 	
		S12[ip] = 0.5*(gradU[1][2][ip] + gradU[2][1][ip]); 
		Smag[ip] = sqrt(2.0*(S00[ip]*S00[ip] + S11[ip]*S11[ip] + S22[ip]*S22[ip]) +
			        4.0*(S01[ip]*S01[ip] + S02[ip]*S02[ip] + S12[ip]*S12[ip])); 	
	    }

	    double *M00_2, *M01_2, *M02_2;
	    double         *M11_2, *M12_2;
	    double 	           *M22_2; 

	    cs->c2d->allocY(M00_2);
	    cs->c2d->allocY(M01_2);
	    cs->c2d->allocY(M02_2);
	    cs->c2d->allocY(M11_2);
	    cs->c2d->allocY(M12_2);
	    cs->c2d->allocY(M22_2);

	    FOR_XYZ_YPEN{
		M00_2[ip] = rho[ip]*Smag[ip]*(S00[ip] - (1.0/3.0)*(S00[p] + S11[ip] + S22[ip]));
		M11_2[ip] = rho[ip]*Smag[ip]*(S11[ip] - (1.0/3.0)*(S00[p] + S11[ip] + S22[ip]));
		M22_2[ip] = rho[ip]*Smag[ip]*(S22[ip] - (1.0/3.0)*(S00[p] + S11[ip] + S22[ip]));

		M01_2[ip] = rho[ip]*Smag[ip]*S01[ip];
		M02_2[ip] = rho[ip]*Smag[ip]*S02[ip];
		M12_2[ip] = rho[ip]*Smag[ip]*S12[ip];
	    }

	    //Filter the Mij_2 fields
	    double *M00_2_hat, *M01_2_hat, *M02_2_hat; 	
	    double 	       *M11_2_hat, *M12_2_hat; 	
	    double 	                   *M22_2_hat; 	

	    cs->c2d->allocY(M00_2_hat);
	    cs->c2d->allocY(M01_2_hat);
	    cs->c2d->allocY(M02_2_hat);
	    cs->c2d->allocY(M11_2_hat);
	    cs->c2d->allocY(M12_2_hat);
	    cs->c2d->allocY(M22_2_hat);

	    filterQuantity(M00_2, M00_2_hat);
	    filterQuantity(M01_2, M01_2_hat);
	    filterQuantity(M02_2, M02_2_hat);
	    filterQuantity(M11_2, M11_2_hat);
	    filterQuantity(M12_2, M12_2_hat);
	    filterQuantity(M22_2, M22_2_hat);

	//Do some averaging of the mu_sgs if we want to...
	doAveraging();
    } 
};



#endif
