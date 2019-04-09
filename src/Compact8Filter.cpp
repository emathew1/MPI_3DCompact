#include "Compact8Filter.hpp"

using namespace std;

void Compact8Filter::multRHSPeriodicFilter(double *phi, double *RHSvec){

    double cc0 = a0_8;
    double cc1 = a1_8/2.0;
    double cc2 = a2_8/2.0;
    double cc3 = a3_8/2.0;
    double cc4 = a4_8/2.0;

    for(int ip = 0; ip < N; ip++){
        if(ip == 0){
            RHSvec[ip] = cc4*phi[N-4]  + cc3*phi[N-3]  +
                         cc2*phi[N-2]  + cc1*phi[N-1]  + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4];
        }else if(ip == 1){
            RHSvec[ip] = cc4*phi[N-3]  + cc3*phi[N-2]  +
                         cc2*phi[N-1]  + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4];
        }else if(ip == 2){
            RHSvec[ip] = cc4*phi[N-2]  + cc3*phi[N-1]  +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4];
        }else if(ip == 3){
            RHSvec[ip] = cc4*phi[N-1]  + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4];
        }else if(ip == N-4){
            RHSvec[ip] = cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[0];
        }else if(ip == N-3){
            RHSvec[ip] = cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[0]    +
                         cc4*phi[1];
        }else if(ip == N-2){
            RHSvec[ip] = cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[0]    + cc3*phi[1]    +
                         cc4*phi[2];   
        }else if(ip == N-1){
            RHSvec[ip] = cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[0]    + cc2*phi[1]    + cc3*phi[2]    +
                         cc4*phi[3];
        }else{           
            RHSvec[ip] = cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4];
        }                
        
    }
}

void Compact8Filter::multRHSDirichletFilter(double *phi, double *RHSvec){

    double aa0 = a0_6;
    double aa1 = a1_6/2.0;
    double aa2 = a2_6/2.0;
    double aa3 = a3_6/2.0;

    double cc0 = a0_8;
    double cc1 = a1_8/2.0;
    double cc2 = a2_8/2.0;
    double cc3 = a3_8/2.0;
    double cc4 = a4_8/2.0;

    RHSvec[0] = phi[0];
    RHSvec[1] = FB2_a*phi[0] + FB2_b*phi[1] + FB2_c*phi[2] + FB2_d*phi[3] + FB2_e*phi[4] +
	        FB2_f*phi[5] + FB2_g*phi[6] + FB2_h*phi[7] + FB2_i*phi[8];
    RHSvec[2] = FB3_a*phi[0] + FB3_b*phi[1] + FB3_c*phi[2] + FB3_d*phi[3] + FB3_e*phi[4] +
	        FB3_f*phi[5] + FB3_g*phi[6] + FB3_h*phi[7] + FB3_i*phi[8];
    RHSvec[3] = FB4_a*phi[0] + FB4_b*phi[1] + FB4_c*phi[2] + FB4_d*phi[3] + FB4_e*phi[4] +
	        FB4_f*phi[5] + FB4_g*phi[6] + FB4_h*phi[7] + FB4_i*phi[8];

//    RHSvec[1] = b00*phi[0] + b01*phi[1] + b02*phi[2] + b03*phi[3] + b04*phi[4];
//    RHSvec[2] = c00*phi[0] + c01*phi[1] + c02*phi[2] + c03*phi[3] + c04*phi[4];
//    RHSvec[3] = aa3*phi[0] + aa2*phi[1] + aa1*phi[2] + aa0*phi[3] + aa1*phi[4] + 
//		aa2*phi[5] + aa3*phi[6];
    
    RHSvec[N-1] = phi[N-1];
//    RHSvec[N-2] = b00*phi[N-1] + b01*phi[N-2] + b02*phi[N-3] + b03*phi[N-4] + b04*phi[N-5];
//    RHSvec[N-3] = c00*phi[N-1] + c01*phi[N-2] + c02*phi[N-3] + c03*phi[N-4] + c04*phi[N-5];
//    RHSvec[N-4] = aa3*phi[N-1] + aa2*phi[N-2] + aa1*phi[N-3] + aa0*phi[N-4] + aa1*phi[N-5] + 
//		  aa2*phi[N-6] + aa3*phi[N-7];
    RHSvec[N-2] = FB2_a*phi[N-1] + FB2_b*phi[N-2] + FB2_c*phi[N-3] + FB2_d*phi[N-4] + FB2_e*phi[N-5] +
	          FB2_f*phi[N-6] + FB2_g*phi[N-7] + FB2_h*phi[N-8] + FB2_i*phi[N-9];
    RHSvec[N-3] = FB3_a*phi[N-1] + FB3_b*phi[N-2] + FB3_c*phi[N-3] + FB3_d*phi[N-4] + FB3_e*phi[N-5] +
	          FB3_f*phi[N-6] + FB3_g*phi[N-7] + FB3_h*phi[N-8] + FB3_i*phi[N-9];
    RHSvec[N-4] = FB4_a*phi[N-1] + FB4_b*phi[N-2] + FB4_c*phi[N-3] + FB4_d*phi[N-4] + FB4_e*phi[N-5] +
	          FB4_f*phi[N-6] + FB4_g*phi[N-7] + FB4_h*phi[N-8] + FB4_i*phi[N-9];


    for(int ip = 4; ip < N-4; ip++){
        RHSvec[ip] = cc4*phi[ip-4] + cc3*phi[ip-3] +
                     cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                     cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                     cc4*phi[ip+4];
    }
}

void Compact8Filter::FilterPeriodic(double *phi, double *phiF){

    double RHSvec[N];
    multRHSPeriodicFilter(phi, RHSvec);
    cyclic(offlowerF, diagF, offupperF, alphaF, alphaF, RHSvec, N, phiF);

}

void Compact8Filter::FilterDirichlet(double *phi, double *phiF){
    

    double RHSvec[N];
    double *work = new double[N];

    multRHSDirichletFilter(phi, RHSvec);
    solveTri(offlowerF, diagF, offupperF, RHSvec, phiF, work, N);

    delete[] work;

}

void Compact8Filter::compactFilter(double *phi, double *phiF){

    if(bcType == Options::PERIODIC_SOLVE){
	FilterPeriodic(phi, phiF);
    }else if(bcType == Options::DIRICHLET_SOLVE){
	FilterDirichlet(phi, phiF);
    }

}

void Compact8Filter::filterField(double *dataIn, double *dataOut){


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
