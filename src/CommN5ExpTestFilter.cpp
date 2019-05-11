#include "CommN5ExpTestFilter.hpp"

using namespace std;

void CommN5ExpTestFilter::FilterPeriodic(double *phi, double *phiF){


    phiF[0] = w3*phi[N-3] + w2*phi[N-2] + w1*phi[N-1] + w0*phi[0] +
	      w1*phi[1]   + w2*phi[2]   + w3*phi[3]; 

    phiF[1] = w3*phi[N-2] + w2*phi[N-1] + w1*phi[0]   + w0*phi[1] +
	      w1*phi[2]   + w2*phi[3]   + w3*phi[4];

    phiF[2] = w3*phi[N-1] + w2*phi[0]   + w1*phi[1]   + w0*phi[2] +
	      w1*phi[3]   + w2*phi[4]   + w3*phi[5];



    for(int ip = 3; ip < N-3; ip++){ 
	phiF[ip] = w3*phi[ip-3] + w2*phi[ip-2] + w1*phi[ip-1] + w0*phi[ip] +
		   w1*phi[ip+1] + w2*phi[ip+2] + w3*phi[ip+3];
    }

    phiF[N-1] = w3*phi[N-4] + w2*phi[N-3] + w1*phi[N-2] + w0*phi[N-1] + 
	       w1*phi[0]   + w2*phi[1]   + w3*phi[2];

    phiF[N-2] = w3*phi[N-5] + w2*phi[N-4] + w1*phi[N-3] + w0*phi[N-2] +
	       w1*phi[N-1] + w2*phi[0]   + w3*phi[1];

    phiF[N-3] = w3*phi[N-6] + w2*phi[N-5] + w1*phi[N-4] + w0*phi[N-3] + 
	       w1*phi[N-2] + w2*phi[N-1] + w3*phi[0];


    cout << "HERE" << endl;
}

void CommN5ExpTestFilter::FilterDirichlet(double *phi, double *phiF){
 
    phiF[0] = w0_n1*phi[0] + w1_n1*phi[1] + w2_n1*phi[2] + w3_n1*phi[3] +
	      w4_n1*phi[4] + w5_n1*phi[5]; 

    phiF[1] = wn1_n2*phi[0] + w0_n2*phi[1] + w1_n2*phi[2] + w2_n2*phi[3] +
	      w3_n2*phi[4]  + w4_n2*phi[5]; 

    phiF[2] = wn2_n3*phi[0] + wn1_n3*phi[1] + w0_n3*phi[2] + w1_n3*phi[3] +
	      w2_n3*phi[4]  + w3_n3*phi[5]; 

    for(int ip = 3; ip < N-3; ip++){ 
	phiF[ip] = w3*phi[ip-3] + w2*phi[ip-2] + w1*phi[ip-1] + w0*phi[ip] +
		   w1*phi[ip+1] + w2*phi[ip+2] + w3*phi[ip+3];
    }

    phiF[N-1] = w0_n1*phi[N-1] + w1_n1*phi[N-2] + w2_n1*phi[N-3] + w3_n1*phi[N-4] +
	        w4_n1*phi[N-5] + w5_n1*phi[N-6]; 

    phiF[N-2] = wn1_n2*phi[N-1] + w0_n2*phi[N-2] + w1_n2*phi[N-3] + w2_n2*phi[N-4] +
	        w3_n2*phi[N-5]  + w4_n2*phi[N-6]; 

    phiF[N-3] = wn2_n3*phi[N-1] + wn1_n3*phi[N-2] + w0_n3*phi[N-3] + w1_n3*phi[N-4] +
	        w2_n3*phi[N-5]  + w3_n3*phi[N-6]; 

}

void CommN5ExpTestFilter::Filter(double *phi, double *phiF){

    if(bcType == Options::PERIODIC_SOLVE){
	FilterPeriodic(phi, phiF);
    }else if(bcType == Options::DIRICHLET_SOLVE){
	FilterDirichlet(phi, phiF);
    }

}

void CommN5ExpTestFilter::filterField(double *dataIn, double *dataOut){


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
		    Filter(dataInLocal, dataOutLocal);
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
                    Filter(dataInLocal, dataOutLocal);
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
		    Filter(dataInLocal, dataOutLocal);
		}
	    }
        }

    }

}
