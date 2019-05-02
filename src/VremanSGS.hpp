#ifndef _CVREMANSGSH_
#define _CVREMANSGSH_

#include "Macros.hpp"
#include "AbstractFilter.hpp"
#include "AbstractDerivatives.hpp"
#include "Compact10Filter.hpp"

class VremanSGS: public AbstractSGS{

  public:

    double c;

    VremanSGS(AbstractCSolver *cs){
	
	this->mpiRank = cs->mpiRank;
	this->cs      = cs;

	Nx = cs->dom->gNx;
	Ny = cs->dom->gNy;
	Nz = cs->dom->gNz;
	N  = cs->dom->gN;
        cs->dom->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);

	musgsAvgType = cs->opt->musgsAvgType;

	//Vreman SGS Coefficient...
	c = 0.07;

	//initialize mu_sgs
	cs->c2d->allocY(mu_sgs);

	//Do some explicit filtering, just need to do something to smooth out the mu_sgs field
	//Could probably do some less complicated filtering or averaging if need be
	filtX = new Compact10Filter(0.0, cs->dom, cs->bc, cs->bc->bcXType, AbstractDerivatives::DIRX);
	filtY = new Compact10Filter(0.0, cs->dom, cs->bc, cs->bc->bcYType, AbstractDerivatives::DIRY);
	filtZ = new Compact10Filter(0.0, cs->dom, cs->bc, cs->bc->bcZType, AbstractDerivatives::DIRZ);
    }
   
    void getSGSViscosity(double *gradU[3][3], double *rho, double *rhoU, double *rhoV, double *rhoW, double *rhoE){

    	//From Vreman (2004)
	FOR_XYZ_YPEN{

	    //There's also an anisotropic form, may need some love to work correctly with curvilinear mesh...
	    double vol = cs->dom->dx*cs->dom->dy*cs->dom->dz/cs->msh->J[ip];
   	    double del = pow(vol, 1.0/3.0);
	    double d1 = del*del;
	    double d2 = d1;
	    double d3 = d1; 
	
	    double d1v1 = gradU[0][0][ip];
	    double d2v1 = gradU[0][1][ip];
	    double d3v1 = gradU[0][2][ip];

	    double d1v2 = gradU[1][0][ip];
	    double d2v2 = gradU[1][1][ip];
	    double d3v2 = gradU[1][2][ip];

	    double d1v3 = gradU[2][0][ip];
	    double d2v3 = gradU[2][1][ip];
	    double d3v3 = gradU[2][2][ip];

	    double b11 = d1*d1v1*d1v1+d2*d2v1*d2v1+d3*d3v1*d3v1;
	    double b12 = d1*d1v1*d1v2+d2*d2v1*d2v2+d3*d3v1*d3v2;
	    double b13 = d1*d1v1*d1v3+d2*d2v1*d2v3+d3*d3v1*d3v3;
	    double b22 = d1*d1v2*d1v2+d2*d2v2*d2v2+d3*d3v2*d3v2;
	    double b23 = d1*d1v2*d1v3+d2*d2v2*d2v3+d3*d3v2*d3v3;
	    double b33 = d1*d1v3*d1v3+d2*d2v3*d2v3+d3*d3v3*d3v3;
	
	    double abeta = d1v1*d1v1 + d1v2*d1v2 + d1v3*d1v3 +
	    		   d2v1*d2v1 + d2v2*d2v2 + d2v3*d2v3 +
			   d3v1*d3v1 + d2v3*d2v3 + d3v3*d3v3;

	    double bbeta = b11*b22-(b12*b12)+b11*b33-(b13*b13)+b22*b33-(b23*b23);

	    double eps = 1e-12;

	    if(bbeta < eps){
		bbeta = eps;
	    }
	    
	    double ss = sqrt(bbeta/(abeta+eps));

	    mu_sgs[ip] = c*rho[ip]*ss;

	    //Should we clip mu_sgs at the high end too?
	    mu_sgs[ip] = min(mu_sgs[ip], 0.05*d1*sqrt(2.0*abeta));
	}

	//Do some averaging of the mu_sgs if we want to...
	doAveraging();
    } 
};



#endif
