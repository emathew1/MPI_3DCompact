#include "CD2.hpp"

void CD2::calc1stDerivField(double *dataIn, double *dataOut){


    if(currentDir == DIRX){

	FOR_Z_XPEN{
	    FOR_Y_XPEN{
		double *dataInLocal, *dataOutLocal;
	        int ii = k*pxSize[0]*pxSize[1] + j*pxSize[0];
		dataInLocal  = &dataIn[ii];
		dataOutLocal = &dataOut[ii];
		calc1stDeriv(dataInLocal, dataOutLocal);
	    }
	}


    }else if(currentDir == DIRY){

	FOR_X_YPEN{
	    FOR_Z_YPEN{
		double *dataInLocal, *dataOutLocal;
	        int ii = i*pySize[2]*pySize[1] + k*pySize[1];
		dataInLocal  = &dataIn[ii];
		dataOutLocal = &dataOut[ii];
		calc1stDeriv(dataInLocal, dataOutLocal);
	    }
	}

    }else if(currentDir == DIRZ){

	FOR_Y_ZPEN{
	    FOR_X_ZPEN{
		double *dataInLocal, *dataOutLocal;
	        int ii = j*pzSize[2]*pzSize[0] + i*pzSize[2];
		dataInLocal  = &dataIn[ii];
		dataOutLocal = &dataOut[ii];
		calc1stDeriv(dataInLocal, dataOutLocal);
	    }
	}

    }

}

//TPB: Transform Periodic Boundaries
void CD2::calc1stDerivField_TPB(double *dataIn, double *dataOut, double *Nm2, double *Nm1, double *Np1, double *Np2){

    if(currentDir == DIRX){

	FOR_Z_XPEN{
	    FOR_Y_XPEN{
	        int ii = k*pxSize[0]*pxSize[1] + j*pxSize[0];

		int iii = k*pxSize[1] + j;

		Nm1val = Nm1[iii];
		Np1val = Np1[iii];

		double *dataInLocal, *dataOutLocal;
		dataInLocal  = &dataIn[ii];
		dataOutLocal = &dataOut[ii];
		calc1stDeriv_TPB(dataInLocal, dataOutLocal);
	    }
	}


    }else if(currentDir == DIRY){

	FOR_X_YPEN{
	    FOR_Z_YPEN{
	        int ii = i*pySize[2]*pySize[1] + k*pySize[1];

		int iii = i*pySize[2] + k;

		Nm1val = Nm1[iii];
		Np1val = Np1[iii];

		double *dataInLocal, *dataOutLocal;
		dataInLocal  = &dataIn[ii];
		dataOutLocal = &dataOut[ii];
		calc1stDeriv_TPB(dataInLocal, dataOutLocal);
	    }
	}

    }else if(currentDir == DIRZ){

	FOR_Y_ZPEN{
	    FOR_X_ZPEN{
	        int ii = j*pzSize[2]*pzSize[0] + i*pzSize[2];

		int iii = j*pzSize[0] + i;

		Nm1val = Nm1[iii];
		Np1val = Np1[iii];

		double *dataInLocal, *dataOutLocal;
		dataInLocal  = &dataIn[ii];
		dataOutLocal = &dataOut[ii];
		calc1stDeriv_TPB(dataInLocal, dataOutLocal);
	    }
	}

    }

}


void CD2::calc1stDeriv(double *phi, double *dphi){

    if(bcType == Options::PERIODIC_SOLVE){
	Compact1stPeriodic(phi, dphi);	
    }else if(bcType == Options::DIRICHLET_SOLVE){
	Compact1stDirichlet(phi, dphi);	
    }

}

void CD2::calc1stDeriv_TPB(double *phi, double *dphi){
    if(bcType == Options::PERIODIC_SOLVE){
	Compact1stPeriodic_TPB(phi, dphi);	
    }else if(bcType == Options::DIRICHLET_SOLVE){
	cout << "CANT USE calc1stDeriv_TPB for non-periodic boundary conditions!!" << endl;
    }
}


void CD2::Compact1stPeriodic(double *phi, double *dphidy){

    dphidy[0]   = (a*phi[N-1] + c*phi[1])/dd;
    dphidy[N-1] = (a*phi[N-2] + c*phi[0])/dd;

    for(int ip = 1; ip < N-1; ip++){
	dphidy[ip] = (a*phi[ip+1] + c*phi[ip-1])/dd;
    } 

    //double RHSvec[N];
    //multRHS1stDerivPeriodic(dd, phi, N, RHSvec);
    //cyclicPenta(offlower2_1D, offlower_1D, diag_1D, offupper_1D, offupper2_1D, RHSvec, cpvec, dphidy, N);

}

void CD2::Compact1stPeriodic_TPB(double *phi, double *dphidy){

    dphidy[0]   = (a*Nm1val   + c*phi[1])/dd;
    dphidy[N-1] = (a*phi[N-2] + c*Np1val)/dd;

    for(int ip = 1; ip < N-1; ip++){
	dphidy[ip] = (a*phi[ip+1] + c*phi[ip-1])/dd;
    } 


/*
    double RHSvec[N];
    multRHS1stDerivPeriodic_TPB(dd, phi, N, RHSvec);
    cyclicPenta(offlower2_1D, offlower_1D, diag_1D, offupper_1D, offupper2_1D, RHSvec, cpvec, dphidy, N);
*/

}

void CD2::Compact1stDirichlet(double *phi, double *dphidy){

    dphidy[0]   = (aa*phi[0]   + bb*phi[1]   + cc*phi[2]   )/dd;
    dphidy[N-1] = (aa*phi[N-1] + bb*phi[N-2] + cc*phi[N-3] )/dd;

    for(int ip = 1; ip < N-1; ip++){
	dphidy[ip] = (a*phi[ip+1] + c*phi[ip-1])/dd;
    } 


/*
    double RHSvec[N];
    double work1[N];
    double work2[N];
    double work3[N];

    multRHS1stDerivDirichlet(dd, phi, N, RHSvec);
    solvePenta(offlower2_1D, offlower_1D, diag_1D, offupper_1D, offupper2_1D, RHSvec, dphidy, work1, work2, work3, N);
*/
}

double CD2::calcNeumann(double *f){
    return (f[0]*4.0 - f[1]*1.0)/3.0;   
}
