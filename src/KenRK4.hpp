#ifndef _CKENRK4H_
#define _CKENRK4H_

#include "Utils.hpp"
#include "AbstractRK.hpp"


//Low Storage, five stage RK4 (Kennedy et al., 1999)

class KenRK4:public AbstractRK{

    public:

	double a1, a2, a3, a4, a5;
	double b1, b2, b3, b4, b5;

        int pxSize[3], pySize[3], pzSize[3]; 
        int pxStart[3], pyStart[3], pzStart[3];
        int pxEnd[3], pyEnd[3], pzEnd[3];

	int baseDirection;
	
	KenRK4(AbstractCSolver *cs){

	    this->cs = cs;
	    baseDirection = cs->baseDirection;
	    mpiRank = cs->mpiRank;
	

	    a1 = 0.0;
	    a2 = -6234157559845.0/12983515589748.0;
	    a3 = -6194124222391.0/4410992767914.0;
	    a4 = -31623096876824.0/15682348800105.0;
	    a5 = -12251185447671.0/11596622555746.0;

	    b1 = 494393426753.0/4806282396855.0;
	    b2 = 4047970641027.0/5463924506627.0;
	    b3 = 9795748752853.0/13190207949281.0;
	    b4 = 4009051133189.0/8539092990294.0;
	    b5 = 1348533437543.0/7166442652324.0;

	    cs->dom->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);
 
	}


        void executeSolverLoop();
        void updateConservedData();

};


void KenRK4::updateConservedData(){
 
    if(cs->useTiming) cs->ft1 = MPI_Wtime();


    if(cs->rkStep == 1){

	if(baseDirection == 1){
            FOR_XYZ_YPEN{
		cs->rhok[ip]  = cs->rho1[ip]  + b1*cs->rhok2[ip]; 
                cs->rhoUk[ip] = cs->rhoU1[ip] + b1*cs->rhoUk2[ip];
                cs->rhoVk[ip] = cs->rhoV1[ip] + b1*cs->rhoVk2[ip];
                cs->rhoWk[ip] = cs->rhoW1[ip] + b1*cs->rhoWk2[ip];
                cs->rhoEk[ip] = cs->rhoE1[ip] + b1*cs->rhoEk2[ip];
	
		cs->rho1[ip]   = cs->rhok2[ip]; 	
		cs->rhoU1[ip]  = cs->rhoUk2[ip]; 	
		cs->rhoV1[ip]  = cs->rhoVk2[ip]; 	
		cs->rhoW1[ip]  = cs->rhoWk2[ip]; 	
		cs->rhoE1[ip]  = cs->rhoEk2[ip]; 	
            }
	}

    }else if(cs->rkStep == 2){
	if(baseDirection == 1){
            FOR_XYZ_YPEN{
                cs->rho1[ip]  = cs->rhok2[ip]  + a2*cs->rho1[ip]; 
                cs->rhoU1[ip] = cs->rhoUk2[ip] + a2*cs->rhoU1[ip];
                cs->rhoV1[ip] = cs->rhoVk2[ip] + a2*cs->rhoV1[ip];
                cs->rhoW1[ip] = cs->rhoWk2[ip] + a2*cs->rhoW1[ip];
                cs->rhoE1[ip] = cs->rhoEk2[ip] + a2*cs->rhoE1[ip];

		cs->rhok[ip]  += b2*cs->rho1[ip]; 
                cs->rhoUk[ip] += b2*cs->rhoU1[ip];
                cs->rhoVk[ip] += b2*cs->rhoV1[ip];
                cs->rhoWk[ip] += b2*cs->rhoW1[ip];
                cs->rhoEk[ip] += b2*cs->rhoE1[ip];
            }
	}

    }else if(cs->rkStep == 3){

	if(baseDirection == 1){
            FOR_XYZ_YPEN{
                cs->rho1[ip]  = cs->rhok2[ip]  + a3*cs->rho1[ip]; 
                cs->rhoU1[ip] = cs->rhoUk2[ip] + a3*cs->rhoU1[ip];
                cs->rhoV1[ip] = cs->rhoVk2[ip] + a3*cs->rhoV1[ip];
                cs->rhoW1[ip] = cs->rhoWk2[ip] + a3*cs->rhoW1[ip];
                cs->rhoE1[ip] = cs->rhoEk2[ip] + a3*cs->rhoE1[ip];

		cs->rhok[ip]  += b3*cs->rho1[ip]; 
                cs->rhoUk[ip] += b3*cs->rhoU1[ip];
                cs->rhoVk[ip] += b3*cs->rhoV1[ip];
                cs->rhoWk[ip] += b3*cs->rhoW1[ip];
                cs->rhoEk[ip] += b3*cs->rhoE1[ip];
	    }
	}

    }else if(cs->rkStep == 4){
	if(baseDirection == 1){
            FOR_XYZ_YPEN{
                cs->rho1[ip]  = cs->rhok2[ip]  + a4*cs->rho1[ip]; 
                cs->rhoU1[ip] = cs->rhoUk2[ip] + a4*cs->rhoU1[ip];
                cs->rhoV1[ip] = cs->rhoVk2[ip] + a4*cs->rhoV1[ip];
                cs->rhoW1[ip] = cs->rhoWk2[ip] + a4*cs->rhoW1[ip];
                cs->rhoE1[ip] = cs->rhoEk2[ip] + a4*cs->rhoE1[ip];

		cs->rhok[ip]  += b4*cs->rho1[ip]; 
                cs->rhoUk[ip] += b4*cs->rhoU1[ip];
                cs->rhoVk[ip] += b4*cs->rhoV1[ip];
                cs->rhoWk[ip] += b4*cs->rhoW1[ip];
                cs->rhoEk[ip] += b4*cs->rhoE1[ip];

	    }
	}
    }else if(cs->rkStep == 5){
	if(baseDirection == 1){
	    FOR_XYZ_YPEN{
                cs->rho1[ip]  = cs->rhok2[ip]  + a5*cs->rho1[ip]; 
                cs->rhoU1[ip] = cs->rhoUk2[ip] + a5*cs->rhoU1[ip];
                cs->rhoV1[ip] = cs->rhoVk2[ip] + a5*cs->rhoV1[ip];
                cs->rhoW1[ip] = cs->rhoWk2[ip] + a5*cs->rhoW1[ip];
                cs->rhoE1[ip] = cs->rhoEk2[ip] + a5*cs->rhoE1[ip];

		cs->rho2[ip]  = cs->rhok[ip]  + b5*cs->rho1[ip]; 
                cs->rhoU2[ip] = cs->rhoUk[ip] + b5*cs->rhoU1[ip];
                cs->rhoV2[ip] = cs->rhoVk[ip] + b5*cs->rhoV1[ip];
                cs->rhoW2[ip] = cs->rhoWk[ip] + b5*cs->rhoW1[ip];
                cs->rhoE2[ip] = cs->rhoEk[ip] + b5*cs->rhoE1[ip];

	    }
	}
    }

    if(cs->useTiming){
        cs->ft2 = MPI_Wtime();
        IF_RANK0 cout << " > updateCons Timing: " << setw(6)  << (int)((cs->ft2-cs->ft1)*1000) << "ms" << endl;
    }

}

void KenRK4::executeSolverLoop(){

    while(cs->endFlag == false){

	cs->rkLast = false;
	
	cs->preStep();

	//Step 1
	cs->rkStep = 1;
	
	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();

	updateConservedData();

	cs->updateData();

	//Step 2
	cs->rkStep = 2;

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();

	//Step 3
	cs->rkStep = 3;

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();


	//Step 4
	cs->rkStep = 4;

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();


	//Step 5
	cs->rkStep = 5;
	cs->rkLast = true;	

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();

	cs->postStep();	

    }
}
        
#endif
