#include "Compact10Filter.hpp"

using namespace std;

void Compact10Filter::multRHSPeriodicFilter(double *phi, double *RHSvec){

    double cc0 = a0_10;
    double cc1 = a1_10/2.0;
    double cc2 = a2_10/2.0;
    double cc3 = a3_10/2.0;
    double cc4 = a4_10/2.0;
    double cc5 = a5_10/2.0;

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

void Compact10Filter::multRHSDirichletFilter(double *phi, double *RHSvec){

    double cc0 = a0_10;
    double cc1 = a1_10/2.0;
    double cc2 = a2_10/2.0;
    double cc3 = a3_10/2.0;
    double cc4 = a4_10/2.0;
    double cc5 = a5_10/2.0;


    RHSvec[0] = FB1_a*phi[0] + FB1_b*phi[1] + FB1_c*phi[2] + FB1_d*phi[3] + FB1_e*phi[4] +
	        FB1_f*phi[5] + FB1_g*phi[6] + FB1_h*phi[7] + FB1_i*phi[8] + FB1_j*phi[9] +
		FB1_k*phi[10];
    RHSvec[1] = FB2_a*phi[0] + FB2_b*phi[1] + FB2_c*phi[2] + FB2_d*phi[3] + FB2_e*phi[4] +
	        FB2_f*phi[5] + FB2_g*phi[6] + FB2_h*phi[7] + FB2_i*phi[8] + FB2_j*phi[9] +
		FB2_k*phi[10];
    RHSvec[2] = FB3_a*phi[0] + FB3_b*phi[1] + FB3_c*phi[2] + FB3_d*phi[3] + FB3_e*phi[4] +
	        FB3_f*phi[5] + FB3_g*phi[6] + FB3_h*phi[7] + FB3_i*phi[8] + FB3_j*phi[9] +
		FB3_k*phi[10];
    RHSvec[3] = FB4_a*phi[0] + FB4_b*phi[1] + FB4_c*phi[2] + FB4_d*phi[3] + FB4_e*phi[4] +
	        FB4_f*phi[5] + FB4_g*phi[6] + FB4_h*phi[7] + FB4_i*phi[8] + FB4_j*phi[9] +
		FB4_k*phi[10];
    RHSvec[4] = FB5_a*phi[0] + FB5_b*phi[1] + FB5_c*phi[2] + FB5_d*phi[3] + FB5_e*phi[4] +
	        FB5_f*phi[5] + FB5_g*phi[6] + FB5_h*phi[7] + FB5_i*phi[8] + FB5_j*phi[9] +
		FB5_k*phi[10];
	
    RHSvec[N-1] = FB1_a*phi[N-1] + FB1_b*phi[N-2] + FB1_c*phi[N-3] + FB1_d*phi[N-4] + FB1_e*phi[N-5] +
	          FB1_f*phi[N-6] + FB1_g*phi[N-7] + FB1_h*phi[N-8] + FB1_i*phi[N-9] + FB1_j*phi[N-10]+
		  FB1_k*phi[N-11];
    RHSvec[N-2] = FB2_a*phi[N-1] + FB2_b*phi[N-2] + FB2_c*phi[N-3] + FB2_d*phi[N-4] + FB2_e*phi[N-5] +
	          FB2_f*phi[N-6] + FB2_g*phi[N-7] + FB2_h*phi[N-8] + FB2_i*phi[N-9] + FB2_j*phi[N-10]+
		  FB2_k*phi[N-11];
    RHSvec[N-3] = FB3_a*phi[N-1] + FB3_b*phi[N-2] + FB3_c*phi[N-3] + FB3_d*phi[N-4] + FB3_e*phi[N-5] +
	          FB3_f*phi[N-6] + FB3_g*phi[N-7] + FB3_h*phi[N-8] + FB3_i*phi[N-9] + FB3_j*phi[N-10]+
		  FB3_k*phi[N-11];
    RHSvec[N-4] = FB4_a*phi[N-1] + FB4_b*phi[N-2] + FB4_c*phi[N-3] + FB4_d*phi[N-4] + FB4_e*phi[N-5] +
	          FB4_f*phi[N-6] + FB4_g*phi[N-7] + FB4_h*phi[N-8] + FB4_i*phi[N-9] + FB4_j*phi[N-10]+
		  FB4_k*phi[N-11];
    RHSvec[N-5] = FB5_a*phi[N-1] + FB5_b*phi[N-2] + FB5_c*phi[N-3] + FB5_d*phi[N-4] + FB5_e*phi[N-5] +
	          FB5_f*phi[N-6] + FB5_g*phi[N-7] + FB5_h*phi[N-8] + FB5_i*phi[N-9] + FB5_j*phi[N-10]+
		  FB5_k*phi[N-11];

    for(int ip = 5; ip < N-5; ip++){
        RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                     cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                     cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                     cc4*phi[ip+4] + cc5*phi[ip+5];
    }
}

void Compact10Filter::FilterPeriodic(double *phi, double *phiF){

    double RHSvec[N];
    multRHSPeriodicFilter(phi, RHSvec);
    cyclic(offlowerF, diagF, offupperF, alphaF, alphaF, RHSvec, N, phiF);

}

void Compact10Filter::FilterDirichlet(double *phi, double *phiF){
    

    double RHSvec[N];
    double *work = new double[N];

    multRHSDirichletFilter(phi, RHSvec);
    solveTri(offlowerF, diagF, offupperF, RHSvec, phiF, work, N);

    delete[] work;

}

void Compact10Filter::compactFilter(double *phi, double *phiF){

    if(bcType == Options::PERIODIC_SOLVE){
	FilterPeriodic(phi, phiF);
    }else if(bcType == Options::DIRICHLET_SOLVE){
	FilterDirichlet(phi, phiF);
    }

}

void Compact10Filter::filterField(double *dataIn, double *dataOut){


    if(currentDir == AbstractDerivatives::DIRX){
	FOR_Z_XPEN{
	    FOR_Y_XPEN{
                double *dataInLocal, *dataOutLocal;
                int ii = k*pxSize[0]*pxSize[1] + j*pxSize[0];
		dataInLocal  = &dataIn[ii];
                dataOutLocal = &dataOut[ii];
		compactFilter(dataInLocal, dataOutLocal);
	    }
	}

    }else if(currentDir == AbstractDerivatives::DIRY){
        FOR_X_YPEN{
            FOR_Z_YPEN{
                double *dataInLocal, *dataOutLocal;
                int ii = i*pySize[2]*pySize[1] + k*pySize[1];
                dataInLocal  = &dataIn[ii];
                dataOutLocal = &dataOut[ii];
                compactFilter(dataInLocal, dataOutLocal);
            }
        }

    }else if(currentDir == AbstractDerivatives::DIRZ){
        FOR_Y_ZPEN{
            FOR_X_ZPEN{
                double *dataInLocal, *dataOutLocal;
                int ii = j*pzSize[2]*pzSize[0] + i*pzSize[2];
                dataInLocal  = &dataIn[ii];
                dataOutLocal = &dataOut[ii];
                compactFilter(dataInLocal, dataOutLocal);
            }
        }

    }

}
