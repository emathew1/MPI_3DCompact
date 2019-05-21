#ifndef _CDSMSGSH_
#define _CDSMSGSH_

#include "Macros.hpp"
#include "AbstractFilter.hpp"
#include "AbstractDerivatives.hpp"
#include "CommN5ExpTestFilter.hpp"
#include "CommN5PadeTestFilter.hpp"

class DSMSGS: public AbstractSGS{

  public:

    double *S00, *S01, *S02;
    double       *S11, *S12;
    double 	       *S22; 
    double *Smag;	    
    double *rho_hat;
    double *Smag_hat;

    double *rhoU_hat, *rhoV_hat, *rhoW_hat;

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
	cs->c2d->allocY(taukk);
	cs->c2d->allocY(Prt_sgs);

	//This should be an input option...
	int filterType = 1;
	double testFilterRatioSquare = 2.0*2.0;

        if(filterType == 1){
	    filtX = new CommN5PadeTestFilter(cs->dom, cs->bc, cs->bc->bcXType, AbstractDerivatives::DIRX);
	    filtY = new CommN5PadeTestFilter(cs->dom, cs->bc, cs->bc->bcYType, AbstractDerivatives::DIRY);
	    filtZ = new CommN5PadeTestFilter(cs->dom, cs->bc, cs->bc->bcZType, AbstractDerivatives::DIRZ);
	}else if(filterType == 2{
	    filtX = new CommN5ExpTestFilter(cs->dom, cs->bc, cs->bc->bcXType, AbstractDerivatives::DIRX);
	    filtY = new CommN5ExpTestFilter(cs->dom, cs->bc, cs->bc->bcYType, AbstractDerivatives::DIRY);
	    filtZ = new CommN5ExpTestFilter(cs->dom, cs->bc, cs->bc->bcZType, AbstractDerivatives::DIRZ);
	}

        cs->c2d->allocY(S00);
        cs->c2d->allocY(S01);
	cs->c2d->allocY(S02);
	cs->c2d->allocY(S11);
	cs->c2d->allocY(S12);
	cs->c2d->allocY(S22);
	cs->c2d->allocY(Smag);
	cs->c2d->allocY(rho_hat);
	cs->c2d->allocY(Smag_hat);	

  	cs->c2d->allocY(rhoU_hat);	   
  	cs->c2d->allocY(rhoV_hat);	   
  	cs->c2d->allocY(rhoW_hat);	   

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
   
    void calcMuSGSTaukk(double *gradU[3][3], double *rho, double *rhoU, double *rhoV, double *rhoW, double *rhoE){


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
		M00_2[ip] = 2.0*rho[ip]*Smag[ip]*(S00[ip] - (1.0/3.0)*(S00[p] + S11[ip] + S22[ip]));
		M11_2[ip] = 2.0*rho[ip]*Smag[ip]*(S11[ip] - (1.0/3.0)*(S00[p] + S11[ip] + S22[ip]));
		M22_2[ip] = 2.0*rho[ip]*Smag[ip]*(S22[ip] - (1.0/3.0)*(S00[p] + S11[ip] + S22[ip]));

		M01_2[ip] = 2.0*rho[ip]*Smag[ip]*S01[ip];
		M02_2[ip] = 2.0*rho[ip]*Smag[ip]*S02[ip];
		M12_2[ip] = 2.0*rho[ip]*Smag[ip]*S12[ip];
	    }

	    //Filter the Mij_2 fields, put them in Mij containers
	    double *M00, *M01, *M02; 	
	    double 	 *M11, *M12; 	
	    double 	       *M22; 	

	    cs->c2d->allocY(M00);
	    cs->c2d->allocY(M01);
	    cs->c2d->allocY(M02);
	    cs->c2d->allocY(M11);
	    cs->c2d->allocY(M12);
	    cs->c2d->allocY(M22);

	    filterQuantity(M00_2, M00);
	    filterQuantity(M01_2, M01);
	    filterQuantity(M02_2, M02);
	    filterQuantity(M11_2, M11);
	    filterQuantity(M12_2, M12);
	    filterQuantity(M22_2, M22);

	    //Filter the M1 field components
	    double *S00_hat, *S01_hat, *S02_hat;
	    double 	     *S11_hat, *S12_hat;
	    double 		       *S22_hat;

	    cs->c2d->allocY(S00_hat);
	    cs->c2d->allocY(S01_hat);
	    cs->c2d->allocY(S02_hat);
	    cs->c2d->allocY(S11_hat);
	    cs->c2d->allocY(S12_hat);
	    cs->c2d->allocY(S22_hat);

	    filterQuantity(rho, rho_hat);
	    filterQuantity(S00, S00_hat);
	    filterQuantity(S01, S01_hat);
	    filterQuantity(S02, S02_hat);
	    filterQuantity(S11, S11_hat);
	    filterQuantity(S12, S12_hat);
	    filterQuantity(S22, S22_hat);
	    filterQuantity(Smag, Smag_hat);

	    //Finish calculating the Mij tensor
	    FOR_XYZ_YPEN{
		M00[ip] -= 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*(S00_hat[ip] - (1.0/3.0)*(S00_hat[ip] + S11_hat[ip] + S22_hat[ip])); 
		M11[ip] -= 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*(S11_hat[ip] - (1.0/3.0)*(S00_hat[ip] + S11_hat[ip] + S22_hat[ip])); 
		M22[ip] -= 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*(S22_hat[ip] - (1.0/3.0)*(S00_hat[ip] + S11_hat[ip] + S22_hat[ip])); 
	        M01[ip] -= 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*S01_hat[ip];
	        M02[ip] -= 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*S02_hat[ip];
	        M12[ip] -= 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*S12_hat[ip];
	    } 

	    //Calculate the stuff for the Leonard stress terms
	    double *L00_t, *L01_t, *L02_t;
	    double 	   *L11_t, *L12_t;
	    double		   *L22_t;
	    
  	    cs->c2d->allocY(L00_t);	   
  	    cs->c2d->allocY(L01_t);	   
  	    cs->c2d->allocY(L02_t);	   
  	    cs->c2d->allocY(L11_t);	   
  	    cs->c2d->allocY(L12_t);	   
  	    cs->c2d->allocY(L22_t);	   
	
	    FOR_XYZ_YPEN{
		L00_t[ip] = rhoU[ip]*rhoU[ip]/rho[ip];
		L01_t[ip] = rhoU[ip]*rhoV[ip]/rho[ip];
		L02_t[ip] = rhoU[ip]*rhoW[ip]/rho[ip];
		L11_t[ip] = rhoV[ip]*rhoV[ip]/rho[ip];
		L12_t[ip] = rhoV[ip]*rhoW[ip]/rho[ip];
		L22_t[ip] = rhoW[ip]*rhoW[ip]/rho[ip];
	    }

	    double *L00, *L01, *L02;
	    double 	 *L11, *L12;
	    double	       *L22;
	     
	    cs->c2d->allocY(L00);	   
  	    cs->c2d->allocY(L01);	   
  	    cs->c2d->allocY(L02);	   
  	    cs->c2d->allocY(L11);	   
  	    cs->c2d->allocY(L12);	   
  	    cs->c2d->allocY(L22);	   
	    
	    filterQuantity(L00_t, L00);
	    filterQuantity(L01_t, L01);
	    filterQuantity(L02_t, L02);
	    filterQuantity(L11_t, L11);
	    filterQuantity(L12_t, L12);
	    filterQuantity(L22_t, L22);


	    filterQuantity(rhoU, rhoU_hat);
	    filterQuantity(rhoV, rhoV_hat);
	    filterQuantity(rhoW, rhoW_hat);

	    FOR_XYZ_YPEN{
		L00[ip] -= rhoU_hat[ip]*rhoU_hat[ip]/rho_hat[ip];
		L01[ip] -= rhoU_hat[ip]*rhoV_hat[ip]/rho_hat[ip];
		L02[ip] -= rhoU_hat[ip]*rhoW_hat[ip]/rho_hat[ip];
		L11[ip] -= rhoV_hat[ip]*rhoV_hat[ip]/rho_hat[ip];
		L12[ip] -= rhoV_hat[ip]*rhoW_hat[ip]/rho_hat[ip];
		L22[ip] -= rhoW_hat[ip]*rhoW_hat[ip]/rho_hat[ip];
	    } 


	    //Get together components for C_I
	    double *CIdenom, *alpha_hat;

	    cs->c2d->allocY(CIdenom);
	    cs->c2d->allocY(alpha_hat);

	    FOR_XYZ_YPEN{
		//going to store the alpha info in beta first...
		CIdenom[ip] = 2.0*rho[ip]*Smag[ip];
	    }

	    //Do the filtering...
	    filterQuantity(CIdenom, alpha_hat);

	    //Then actually calculate the denomenator 
	    FOR_XYZ_YPEN{
		CIdenom[ip] = 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*Smag_hat[ip];
	    }

	    //Need to do some quantity averaging now in whatever kind of averaging we have planned...
	    double *Lkk, *Lkk_avg, *CIdenom_avg;
	    double *LijMij, *LijMij_avg, *MijMij, *MijMij_avg;
	    cs->c2d->allocY(Lkk); 
	    cs->c2d->allocY(CIdenom_avg);
	    cs->c2d->allocY(LijMij); 
	    cs->c2d->allocY(MijMij); 

	    FOR_XYZ_YPEN{
		Lkk[ip] = L00[ip] + L11[ip] + L22[ip];

		LijMij[ip] = (L00[ip] - (1.0/3.0)*Lkk[ip])*M00[ip] +
			     (L11[ip] - (1.0/3.0)*Lkk[ip])*M11[ip] +
			     (L22[ip] - (1.0/3.0)*Lkk[ip])*M22[ip] +
			      2.0*L01[ip]*M01[ip] +
			      2.0*L02[ip]*M02[ip] +
			      2.0*L12[ip]*M12[ip];
		MijMij[ip] = M00[ip]*M00[ip] + M11[ip]*M11[ip] + M22[ip]*M22[ip] +
			     2.0*M01[ip]*M01[ip] +
			     2.0*M02[ip]*M02[ip] + 
			     2.0*M12[ip]*M12[ip]; 
	    }

	    cs->c2d->allocY(Lkk_avg);
	    doAveraging(Lkk, Lkk_avg);

	    cs->c2d->allocY(CIdenom_avg);
	    doAveraging(CIdenom, CIdenom_avg);

	    cs->c2d->allocY(LijMij_avg);
	    doAveraging(LijMij, LijMij_avg);

	    cs->c2d->allocY(MijMij_avg);
	    doAveraging(MijMij_avg);

	    FOR_XYZ_YPEN{
		mu_sgs[ip] = rho[ip]*Smag[ip]*LijMij_avg[ip]/MijMij_avg[ip];
		taukk[ip]  = rho[ip]*Smag[ip]*Smag[ip]*Lkk_avg[ip]/CIdenom_avg[ip];
	    }

    }

    void calcPrtSGS(double *gradT[3], double *rho, double *rhoU, double *rhoV, double *rhoW, double *T){

    } 
};



#endif
