#ifndef _CDSMSGSH_
#define _CDSMSGSH_

#include "Macros.hpp"
#include "AbstractSGS.hpp"
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

    double *work1,  *work2,  *work3,  *work4,  *work5, *work6;
    double *work7,  *work8,  *work9,  *work10, *work11, *work12;
    double *work13, *work14, *work15, *work16, *work17, *work18;

    double *rhoU_hat, *rhoV_hat, *rhoW_hat;

    double testFilterRatioSquare;
    bool useTaukk, dumpDSMCoeff;
    int filterType;

    DSMSGS(AbstractCSolver *cs, int filterType, bool useTaukk, bool dumpDSMCoeff){
	
	this->mpiRank = cs->mpiRank;
	this->cs      = cs;
	
	this->filterType = filterType;
        this->useTaukk = useTaukk;
	this->dumpDSMCoeff = dumpDSMCoeff;

	Nx = cs->dom->gNx;
	Ny = cs->dom->gNy;
	Nz = cs->dom->gNz;
	N  = cs->dom->gN;
        cs->dom->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);

	musgsAvgType = cs->opt->musgsAvgType;

	//initialize mu_sgs
	cs->c2d->allocY(mu_sgs);
	cs->c2d->allocY(taukk);
	cs->c2d->allocY(k_sgs);

	//This should be an input option...
	testFilterRatioSquare = 2.0*2.0;

        if(filterType == 1){
	    filtX = new CommN5PadeTestFilter(cs->dom, cs->bc, cs->bc->bcXType, AbstractDerivatives::DIRX);
	    filtY = new CommN5PadeTestFilter(cs->dom, cs->bc, cs->bc->bcYType, AbstractDerivatives::DIRY);
	    filtZ = new CommN5PadeTestFilter(cs->dom, cs->bc, cs->bc->bcZType, AbstractDerivatives::DIRZ);
	}else if(filterType == 2){
	    filtX = new CommN5ExpTestFilter(cs->dom, cs->bc, cs->bc->bcXType, AbstractDerivatives::DIRX);
	    filtY = new CommN5ExpTestFilter(cs->dom, cs->bc, cs->bc->bcYType, AbstractDerivatives::DIRY);
	    filtZ = new CommN5ExpTestFilter(cs->dom, cs->bc, cs->bc->bcZType, AbstractDerivatives::DIRZ);
	}else{
	    cout << "ABORT, UNKNOWN TEST FILTER TYPE" << endl;
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

  	cs->c2d->allocY(work1);	   
  	cs->c2d->allocY(work2);	   
  	cs->c2d->allocY(work3);	   
  	cs->c2d->allocY(work4);	   
  	cs->c2d->allocY(work5);	   
  	cs->c2d->allocY(work6);
	   
  	cs->c2d->allocY(work7);	   
  	cs->c2d->allocY(work8);	   
  	cs->c2d->allocY(work9);	   
  	cs->c2d->allocY(work10);	   
  	cs->c2d->allocY(work11);	   
  	cs->c2d->allocY(work12);	   
  	
	cs->c2d->allocY(work13);	   
	cs->c2d->allocY(work14);	   
	cs->c2d->allocY(work15);	   
	cs->c2d->allocY(work16);	   
	cs->c2d->allocY(work17);	   
	cs->c2d->allocY(work18);	   

	//If we want to dump smag coeff ranges...
	if(dumpDSMCoeff){
	    cs->c2d->allocY(C);
	    cs->c2d->allocY(CI);
	    cs->c2d->allocY(Prt);
	}
 
   }

    void filterQuantity(double *phi, double *phiF){

	filtY->filterField(phi, cs->tempY1);
	cs->c2d->transposeY2Z_MajorIndex(cs->tempY1, cs->tempZ1);
	filtZ->filterField(cs->tempZ1, cs->tempZ2);
	cs->c2d->transposeZ2Y_MajorIndex(cs->tempZ2, cs->tempY1);
	cs->c2d->transposeY2X_MajorIndex(cs->tempY1, cs->tempX1);
	filtX->filterField(cs->tempX1, cs->tempX2);
	cs->c2d->transposeX2Y_MajorIndex(cs->tempX2, phiF);
    }
   
    void calcMuSGSTaukk(double *gradU[3][3], double *rho, double *rhoU, double *rhoV, double *rhoW, double *rhoE){


	    double *M00_2, *M01_2, *M02_2;
	    double         *M11_2, *M12_2;
	    double 	           *M22_2; 
	
	    M00_2 = work1;
	    M01_2 = work2;
	    M02_2 = work3;
	    M11_2 = work4;
	    M12_2 = work5;
	    M22_2 = work6;


	    FOR_XYZ_YPEN{
		S00[ip] = gradU[0][0][ip];
		S11[ip] = gradU[1][1][ip];
		S22[ip] = gradU[2][2][ip];
		S01[ip] = 0.5*(gradU[0][1][ip] + gradU[1][0][ip]); 	
		S02[ip] = 0.5*(gradU[0][2][ip] + gradU[2][0][ip]); 	
		S12[ip] = 0.5*(gradU[1][2][ip] + gradU[2][1][ip]); 
		Smag[ip] = sqrt(2.0*(S00[ip]*S00[ip] + S11[ip]*S11[ip] + S22[ip]*S22[ip]) +
			        4.0*(S01[ip]*S01[ip] + S02[ip]*S02[ip] + S12[ip]*S12[ip])); 	

		M00_2[ip] = 2.0*rho[ip]*Smag[ip]*(S00[ip] - (1.0/3.0)*(S00[ip] + S11[ip] + S22[ip]));
		M11_2[ip] = 2.0*rho[ip]*Smag[ip]*(S11[ip] - (1.0/3.0)*(S00[ip] + S11[ip] + S22[ip]));
		M22_2[ip] = 2.0*rho[ip]*Smag[ip]*(S22[ip] - (1.0/3.0)*(S00[ip] + S11[ip] + S22[ip]));

		M01_2[ip] = 2.0*rho[ip]*Smag[ip]*S01[ip];
		M02_2[ip] = 2.0*rho[ip]*Smag[ip]*S02[ip];
		M12_2[ip] = 2.0*rho[ip]*Smag[ip]*S12[ip];
	    }

	    //Filter the Mij_2 fields, put them in Mij containers
	    double *M00, *M01, *M02; 	
	    double 	 *M11, *M12; 	
	    double 	       *M22; 	

	    M00 = work7;
	    M01 = work8;
	    M02 = work9;
	    M11 = work10;
	    M12 = work11;
	    M22 = work12;

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

	    S00_hat = work1;
	    S01_hat = work2;
	    S02_hat = work3;
	    S11_hat = work4;
	    S12_hat = work5;
	    S22_hat = work6;

	    filterQuantity(rho, rho_hat);
	    filterQuantity(S00, S00_hat);
	    filterQuantity(S01, S01_hat);
	    filterQuantity(S02, S02_hat);
	    filterQuantity(S11, S11_hat);
	    filterQuantity(S12, S12_hat);
	    filterQuantity(S22, S22_hat);
	    filterQuantity(Smag, Smag_hat);

	    //Calculate the stuff for the Leonard stress terms
	    double *L00_t, *L01_t, *L02_t;
	    double 	   *L11_t, *L12_t;
	    double		   *L22_t;

	    L00_t = work13;	    
	    L01_t = work14;	    
	    L02_t = work15;	    
	    L11_t = work16;	    
	    L12_t = work17;	    
	    L22_t = work18;	    
	

	    //Finish calculating the Mij tensor
	    FOR_XYZ_YPEN{
		M00[ip] -= 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*(S00_hat[ip] - (1.0/3.0)*(S00_hat[ip] + S11_hat[ip] + S22_hat[ip])); 
		M11[ip] -= 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*(S11_hat[ip] - (1.0/3.0)*(S00_hat[ip] + S11_hat[ip] + S22_hat[ip])); 
		M22[ip] -= 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*(S22_hat[ip] - (1.0/3.0)*(S00_hat[ip] + S11_hat[ip] + S22_hat[ip])); 
	        M01[ip] -= 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*S01_hat[ip];
	        M02[ip] -= 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*S02_hat[ip];
	        M12[ip] -= 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*S12_hat[ip];

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
	     
	    L00 = work1;	   
  	    L01 = work2;	   
  	    L02 = work3;	   
  	    L11 = work4;	   
  	    L12 = work5;	   
  	    L22 = work6;	   
	    
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

	    double *CIdenom, *alpha_hat;
	    if(useTaukk){
	        //Get together components for C_I

	        CIdenom   = work13;
	        alpha_hat = work14;

	        FOR_XYZ_YPEN{
		    //going to store the alpha info in beta first...
		    CIdenom[ip] = 2.0*rho[ip]*Smag[ip]*Smag[ip];
	        }

	        //Do the filtering...
	        filterQuantity(CIdenom, alpha_hat);

	        //Then actually calculate the denomenator 
	        FOR_XYZ_YPEN{
		    CIdenom[ip] = 2.0*testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*Smag_hat[ip] - alpha_hat[ip];
	        }
	    }

	    //Need to do some quantity averaging now in whatever kind of averaging we have planned...
	    double *Lkk, *Lkk_avg, *CIdenom_avg;
	    double *LijMij, *LijMij_avg, *MijMij, *MijMij_avg;
	    Lkk = work15; 
	    CIdenom_avg = work16;
	    LijMij = work17; 
	    MijMij = work18; 

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

	    if(useTaukk){
	        Lkk_avg = work7;
	        doAveraging(Lkk, Lkk_avg);

	        CIdenom_avg = work8;
	        doAveraging(CIdenom, CIdenom_avg);
	    }

	    LijMij_avg = work9;
	    doAveraging(LijMij, LijMij_avg);

	    MijMij_avg = work10;
	    doAveraging(MijMij, MijMij_avg);

	    FOR_XYZ_YPEN{
		mu_sgs[ip] = rho[ip]*Smag[ip]*LijMij_avg[ip]/MijMij_avg[ip];

		if(useTaukk){
		    taukk[ip]  = 2.0*rho[ip]*Smag[ip]*Smag[ip]*Lkk_avg[ip]/CIdenom_avg[ip];
		}else{
		    taukk[ip]  = 0.0;
		}
	    }

	    //do something here to calculate the actual C & CI and dump it out...

	    if(dumpDSMCoeff){	
		FOR_XYZ_YPEN{
	            double vol = cs->dom->dx*cs->dom->dy*cs->dom->dz/cs->msh->J[ip];
   	            double del = pow(vol, 1.0/3.0);
		    C[ip]  = LijMij_avg[ip]/(del*del*MijMij_avg[ip]);
		    if(useTaukk){
		         CI[ip] = Lkk_avg[ip]/(del*del*CIdenom_avg[ip]);
		    }else{
			 CI[ip] = 0.0;
		    }
		}
	    }

    }

    void calcKSGS(double *gradT[3], double *rho, double *rhoU, double *rhoV, double *rhoW, double *T){

	double *N_0, *N_1, *N_2;	
	double *NRHS_0, *NRHS_1, *NRHS_2;	

  	N_0 = work1;	   
  	N_1 = work2;	   
  	N_2 = work3;	   
  	NRHS_0  = work4;	   
  	NRHS_1  = work5;	   
  	NRHS_2  = work6;	   

	double *K_0, *K_1, *K_2;	
	double *KRHS_0, *KRHS_1, *KRHS_2;	

	K_0 = work7;
	K_1 = work8;
	K_2 = work9;
	KRHS_0 = work10; 
	KRHS_1 = work11; 
	KRHS_2 = work12; 

	double *gradT0_hat, *gradT1_hat, *gradT2_hat; 
  	gradT0_hat = work13;	   
  	gradT1_hat = work14;	   
  	gradT2_hat = work15;	   
	double *rhoT, *rhoT_hat;
  	rhoT = work16;	   
  	rhoT_hat = work17;	   

	FOR_XYZ_YPEN{
	    //We'll store the NRHS stuff in N_i until we filter...
	    N_0[ip] = rho[ip]*Smag[ip]*gradT[0][ip];
	    N_1[ip] = rho[ip]*Smag[ip]*gradT[1][ip];
	    N_2[ip] = rho[ip]*Smag[ip]*gradT[2][ip];

	    K_0[ip] = rhoU[ip]*T[ip];
	    K_1[ip] = rhoV[ip]*T[ip];
	    K_2[ip] = rhoW[ip]*T[ip];
	    rhoT[ip] = rho[ip]*T[ip];
	}

	filterQuantity(N_0, NRHS_0);		
	filterQuantity(N_1, NRHS_1);		
	filterQuantity(N_2, NRHS_2);		

	filterQuantity(K_0, KRHS_0);
	filterQuantity(K_1, KRHS_1);
	filterQuantity(K_2, KRHS_2);

	filterQuantity(rhoT, rhoT_hat);

	filterQuantity(gradT[0], gradT0_hat);
	filterQuantity(gradT[1], gradT1_hat);
	filterQuantity(gradT[2], gradT2_hat);

	FOR_XYZ_YPEN{
	    N_0[ip] = testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*gradT0_hat[ip] - NRHS_0[ip]; 
	    N_1[ip] = testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*gradT1_hat[ip] - NRHS_1[ip]; 
	    N_2[ip] = testFilterRatioSquare*rho_hat[ip]*Smag_hat[ip]*gradT2_hat[ip] - NRHS_2[ip]; 

	    K_0[ip] = -rhoU_hat[ip]*rhoT_hat[ip]/rho[ip] + KRHS_0[ip];
	    K_1[ip] = -rhoV_hat[ip]*rhoT_hat[ip]/rho[ip] + KRHS_1[ip];
	    K_2[ip] = -rhoW_hat[ip]*rhoT_hat[ip]/rho[ip] + KRHS_2[ip];
	}

	double *NiNi, *NiKi;
  	NiNi = work4;	   
  	NiKi = work5;	   

	FOR_XYZ_YPEN{
	    NiNi[ip] = N_0[ip]*N_0[ip] + N_1[ip]*N_1[ip] + N_2[ip]*N_2[ip];
	    NiKi[ip] = K_0[ip]*N_0[ip] + K_1[ip]*N_1[ip] + K_2[ip]*N_2[ip];
	}

	double *NiNi_avg, *NiKi_avg;
	NiNi_avg = work6;
	NiKi_avg = work10;

	doAveraging(NiNi, NiNi_avg);
	doAveraging(NiKi, NiKi_avg);

	//Calculate -C/Pr_t = KjNj/NiNi
	FOR_XYZ_YPEN{
	    double C_over_Prt = NiKi_avg[ip]/NiNi_avg[ip];
	    //Should this be negative? Don't think so...
	    k_sgs[ip] = -cs->ig->cp*C_over_Prt*rho[ip]*Smag[ip];
	}

	
	if(dumpDSMCoeff){	
	    FOR_XYZ_YPEN{
	        double vol = cs->dom->dx*cs->dom->dy*cs->dom->dz/cs->msh->J[ip];
   	        double del = pow(vol, 1.0/3.0);

		Prt[ip] = -C[ip]*del*del*NiNi_avg[ip]/NiKi_avg[ip];
	    }
	}

    } 

    void getSGSViscosity(double *gradU[3][3], double *rho, double *rhoU, double *rhoV, double *rhoW, double *rhoE){};
};



#endif
