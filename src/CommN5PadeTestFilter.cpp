#include "CommN5PadeTestFilter.hpp"

using namespace std;

void CommN5PadeTestFilter::multRHSPeriodicFilter(double *phi, double *RHSvec){

    double cc5 = wi5;
    double cc4 = wi4;
    double cc3 = wi3;
    double cc2 = wi2;
    double cc1 = wi1;
    double cc0 = wi0;

    for(int ip = 0; ip < N; ip++){
        if(ip == 0){
            RHSvec[ip] = cc5*phi[N-5]  + cc4*phi[N-4]  + cc3*phi[N-3]  +
                         cc2*phi[N-2]  + cc1*phi[N-1]  + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[ip+5];
        }else if(ip == 1){
            RHSvec[ip] = cc5*phi[N-4]  + cc4*phi[N-3]  + cc3*phi[N-2]  +
                         cc2*phi[N-1]  + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[ip+5];
        }else if(ip == 2){
            RHSvec[ip] = cc5*phi[N-3]  + cc4*phi[N-2]  + cc3*phi[N-1]  +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[ip+5];
        }else if(ip == 3){
            RHSvec[ip] = cc5*phi[N-2]  + cc4*phi[N-1]  + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[ip+5];
        }else if(ip == 4){
            RHSvec[ip] = cc5*phi[N-1]  + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[ip+5];

        }else if(ip == N-5){
            RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[0];
        }else if(ip == N-4){
            RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[0]    + cc5*phi[1];
        }else if(ip == N-3){
            RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[0]    +
                         cc4*phi[1]    + cc5*phi[2];
        }else if(ip == N-2){
            RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[0]    + cc3*phi[1]    +
                         cc4*phi[2]    + cc5*phi[3];   
        }else if(ip == N-1){
            RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[0]    + cc2*phi[1]    + cc3*phi[2]    +
                         cc4*phi[3]    + cc5*phi[4];
        }else{           
            RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[ip+5];
        }                
        
    }
}

void CommN5PadeTestFilter::multRHSDirichletFilter(double *phi, double *RHSvec){

    double cc5 = wi5;
    double cc4 = wi4;
    double cc3 = wi3;
    double cc2 = wi2;
    double cc1 = wi1;
    double cc0 = wi0;

    RHSvec[0] = w0_n1*phi[0] + w1_n1*phi[1] + w2_n1*phi[2] + w3_n1*phi[3] +
	        w4_n1*phi[4] + w5_n1*phi[5]; 

    RHSvec[1] = wn1_n2*phi[0] + w0_n2*phi[1] + w1_n2*phi[2] + w2_n2*phi[3] +
	        w3_n2*phi[4]  + w4_n2*phi[5]; 

    RHSvec[2] = wn2_n3*phi[0] + wn1_n3*phi[1] + w0_n3*phi[2] + w1_n3*phi[3] +
	        w2_n3*phi[4]  + w3_n3*phi[5]; 

    RHSvec[3] = w3*phi[0] + w2*phi[1] + w1*phi[2] + w0*phi[3] + w1*phi[4] +
	        w2*phi[5] + w3*phi[6];

    RHSvec[4] = w3*phi[1] + w2*phi[2] + w1*phi[3] + w0*phi[4] + w1*phi[5] +
	        w2*phi[6] + w3*phi[7];

    for(int ip = 5; ip < N-5; ip++){
        RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                     cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                     cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                     cc4*phi[ip+4] + cc5*phi[ip+5];
    }

    RHSvec[N-1] = w0_n1*phi[N-1] + w1_n1*phi[N-2] + w2_n1*phi[N-3] + w3_n1*phi[N-4] +
	          w4_n1*phi[N-5] + w5_n1*phi[N-6]; 

    RHSvec[N-2] = wn1_n2*phi[N-1] + w0_n2*phi[N-2] + w1_n2*phi[N-3] + w2_n2*phi[N-4] +
	          w3_n2*phi[N-5]  + w4_n2*phi[N-6]; 

    RHSvec[N-3] = wn2_n3*phi[N-1] + wn1_n3*phi[N-2] + w0_n3*phi[N-3] + w1_n3*phi[N-4] +
	          w2_n3*phi[N-5]  + w3_n3*phi[N-6]; 

    RHSvec[N-4] = w3*phi[N-1] + w2*phi[N-2] + w1*phi[N-3] + w0*phi[N-4] +
		  w1*phi[N-5] + w2*phi[N-6] + w3*phi[N-7];
 
    RHSvec[N-5] = w3*phi[N-2] + w2*phi[N-3] + w1*phi[N-4] + w0*phi[N-5] +
		  w1*phi[N-6] + w2*phi[N-7] + w3*phi[N-8];
 

}

void CommN5PadeTestFilter::FilterPeriodic(double *phi, double *phiF){

    double RHSvec[N];
    multRHSPeriodicFilter(phi, RHSvec);
    cyclicPenta(offlowerF2, offlowerF, diagF, offupperF, offupperF2, RHSvec,  cpvec, phiF, N);

}

void CommN5PadeTestFilter::FilterDirichlet(double *phi, double *phiF){
    

    double RHSvec[N];
    double work1[N];
    double work2[N];
    double work3[N];

    multRHSDirichletFilter(phi, RHSvec);
    solvePenta(offlowerF2, offlowerF, diagF, offupperF, offupperF2, RHSvec, phiF, work1, work2, work3, N);

}

void CommN5PadeTestFilter::compactFilter(double *phi, double *phiF){

    if(bcType == Options::PERIODIC_SOLVE){
	FilterPeriodic(phi, phiF);
    }else if(bcType == Options::DIRICHLET_SOLVE){
	FilterDirichlet(phi, phiF);
    }

}

void CommN5PadeTestFilter::filterField(double *dataIn, double *dataOut){


    if(currentDir == AbstractDerivatives::DIRX){
	FOR_Z_XPEN{
    	    FOR_Y_XPEN{
	   	double kGlob = GETGLOBALZIND_XPEN;
	  	double jGlob = GETGLOBALYIND_XPEN;

		bool passFlag = false;
		if(bc->bcZType == Options::DIRICHLET_SOLVE && (kGlob == 0 || kGlob == Nz-1)){
		    passFlag =true;
		}

		if(bc->bcYType == Options::DIRICHLET_SOLVE && (jGlob == 0 || jGlob == Ny-1)){
		    passFlag = true;
		}

    		double *dataInLocal, *dataOutLocal;
                int ii = k*pxSize[0]*pxSize[1] + j*pxSize[0];
	        dataInLocal  = &dataIn[ii];
                dataOutLocal = &dataOut[ii];
		if(passFlag){
	            memcpy(dataOutLocal, dataInLocal, sizeof(double)*pxSize[0]);
		}else{
		    compactFilter(dataInLocal, dataOutLocal);
		}

   	    }
	}

    }else if(currentDir == AbstractDerivatives::DIRY){
        FOR_X_YPEN{
            FOR_Z_YPEN{
	        double iGlob = GETGLOBALXIND_YPEN;
	        double kGlob = GETGLOBALZIND_YPEN;

		bool passFlag = false;
		if(bc->bcXType == Options::DIRICHLET_SOLVE && (iGlob == 0 || iGlob == Nx-1))
		    passFlag = true;

		if(bc->bcZType == Options::DIRICHLET_SOLVE && (kGlob == 0 || kGlob == Nz-1))
		    passFlag = true;


		double *dataInLocal, *dataOutLocal;
                int ii = i*pySize[2]*pySize[1] + k*pySize[1];
                dataInLocal  = &dataIn[ii];
                dataOutLocal = &dataOut[ii];
		if(passFlag){
		    memcpy(dataOutLocal, dataInLocal, sizeof(double)*pySize[1]);
		}else{
                    compactFilter(dataInLocal, dataOutLocal);
		}

	    }
        }

    }else if(currentDir == AbstractDerivatives::DIRZ){
        FOR_Y_ZPEN{
            FOR_X_ZPEN{
	        double jGlob = GETGLOBALYIND_ZPEN;
	        double iGlob = GETGLOBALXIND_ZPEN;

	        bool passFlag = false;	
	     	if(bc->bcYType == Options::DIRICHLET_SOLVE && (jGlob == 0 || jGlob == Ny-1))
		    passFlag = true;

		if(bc->bcXType == Options::DIRICHLET_SOLVE && (iGlob == 0 || iGlob == Nx-1))
		    passFlag = true;

    		double *dataInLocal, *dataOutLocal;    
		int ii = j*pzSize[2]*pzSize[0] + i*pzSize[2];
                dataInLocal  = &dataIn[ii];
                dataOutLocal = &dataOut[ii];

		if(passFlag){
		    memcpy(dataOutLocal, dataInLocal, sizeof(double)*pzSize[2]);
		}else{
		    compactFilter(dataInLocal, dataOutLocal);
		}
	    }
        }

    }

}
