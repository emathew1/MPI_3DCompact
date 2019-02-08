#include "Penta10.hpp"

void Penta10::calc1stDerivField(double *dataIn, double *dataOut){


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
void Penta10::calc1stDerivField_TPB(double *dataIn, double *dataOut, double *Nm3, double *Nm2, double *Nm1, double *Np1, double *Np2, double *Np3){

    if(currentDir == DIRX){

	FOR_Z_XPEN{
	    FOR_Y_XPEN{
	        int ii = k*pxSize[0]*pxSize[1] + j*pxSize[0];

		int iii = k*pxSize[1] + j;

		Nm3val = Nm3[iii];
		Nm2val = Nm2[iii];
		Nm1val = Nm1[iii];
		Np1val = Np1[iii];
		Np2val = Np2[iii];
		Np3val = Np3[iii];

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

		Nm3val = Nm3[iii];
		Nm2val = Nm2[iii];
		Nm1val = Nm1[iii];
		Np1val = Np1[iii];
		Np2val = Np2[iii];
		Np3val = Np3[iii];

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

		Nm3val = Nm3[iii];
		Nm2val = Nm2[iii];
		Nm1val = Nm1[iii];
		Np1val = Np1[iii];
		Np2val = Np2[iii];
		Np3val = Np3[iii];

		double *dataInLocal, *dataOutLocal;
		dataInLocal  = &dataIn[ii];
		dataOutLocal = &dataOut[ii];
		calc1stDeriv_TPB(dataInLocal, dataOutLocal);
	    }
	}

    }

}


void Penta10::calc1stDeriv(double *phi, double *dphi){

    if(bcType == Options::PERIODIC_SOLVE){
	Compact1stPeriodic(phi, dphi);	
    }else if(bcType == Options::DIRICHLET_SOLVE){
	Compact1stDirichlet(phi, dphi);	
    }

}

void Penta10::calc1stDeriv_TPB(double *phi, double *dphi){
    if(bcType == Options::PERIODIC_SOLVE){
	Compact1stPeriodic_TPB(phi, dphi);	
    }else if(bcType == Options::DIRICHLET_SOLVE){
	cout << "CANT USE calc1stDeriv_TPB for non-periodic boundary conditions!!" << endl;
    }
}


void Penta10::multRHS1stDerivPeriodic(double dh, double *phi, int N, double *RHSvec){

    double cc1 = -c/6.0;
    double cc2 = -b/4.0;
    double cc3 = -a/2.0;
    double cc4 =  a/2.0;
    double cc5 =  b/4.0;
    double cc6 =  c/6.0;

    RHSvec[0] = cc1*phi[N-3] + cc2*phi[N-2] + cc3*phi[N-1] + \
                        cc4*phi[1] + cc5*phi[2] + cc6*phi[3];
    RHSvec[1] = cc1*phi[N-2] + cc2*phi[N-1] + cc3*phi[0] + \
                        cc4*phi[2] + cc5*phi[3] + cc6*phi[4];
    RHSvec[2] = cc1*phi[N-1] + cc2*phi[0] + cc3*phi[1] + \
                        cc4*phi[3] + cc5*phi[4] + cc6*phi[5];

    for(int ip = 3; ip < N-3; ip++){
        RHSvec[ip] = cc1*phi[ip-3] + cc2*phi[ip-2] + cc3*phi[ip-1] + \
                        cc4*phi[ip+1] + cc5*phi[ip+2] + cc6*phi[ip+3];
    }

    RHSvec[N-3] = cc1*phi[N-6] + cc2*phi[N-5] + cc3*phi[N-4] + \
                        cc4*phi[N-2] + cc5*phi[N-1] + cc6*phi[0];
    RHSvec[N-2] = cc1*phi[N-5] + cc2*phi[N-4] + cc3*phi[N-3] + \
                        cc4*phi[N-1] + cc5*phi[0] + cc6*phi[1];
    RHSvec[N-1] = cc1*phi[N-4] + cc2*phi[N-3] + cc3*phi[N-2] + \
                        cc4*phi[0] + cc5*phi[1] + cc6*phi[2];

    for(int ip = 0; ip < N; ip++){
        RHSvec[ip] /= dh;
    }

}

void Penta10::multRHS1stDerivPeriodic_TPB(double dh, double *phi, int N, double *RHSvec){
 
    double cc1 = -c/6.0;
    double cc2 = -b/4.0;
    double cc3 = -a/2.0;
    double cc4 =  a/2.0;
    double cc5 =  b/4.0;
    double cc6 =  c/6.0;

    RHSvec[0] = cc1*Nm3val + cc2*Nm2val + cc3*Nm1val + \
                        cc4*phi[1] + cc5*phi[2] + cc6*phi[3];
    RHSvec[1] = cc1*Nm2val + cc2*Nm1val + cc3*phi[0] + \
                        cc4*phi[2] + cc5*phi[3] + cc6*phi[4];
    RHSvec[2] = cc1*Nm1val + cc2*phi[0] + cc3*phi[1] + \
                        cc4*phi[3] + cc5*phi[4] + cc6*phi[5];

    for(int ip = 3; ip < N-3; ip++){
        RHSvec[ip] = cc1*phi[ip-3] + cc2*phi[ip-2] + cc3*phi[ip-1] + \
                        cc4*phi[ip+1] + cc5*phi[ip+2] + cc6*phi[ip+3];

    }

    RHSvec[N-3] = cc1*phi[N-6] + cc2*phi[N-5] + cc3*phi[N-4] + \
                        cc4*phi[N-2] + cc5*phi[N-1] + cc6*Np1val;
    RHSvec[N-2] = cc1*phi[N-5] + cc2*phi[N-4] + cc3*phi[N-3] + \
                        cc4*phi[N-1] + cc5*Np1val + cc6*Np2val;
    RHSvec[N-1] = cc1*phi[N-4] + cc2*phi[N-3] + cc3*phi[N-2] + \
                        cc4*Np1val + cc5*Np2val + cc6*Np3val;

    for(int ip = 0; ip < N; ip++){
        RHSvec[ip] /= dh;
    }

}


void Penta10::multRHS1stDerivDirichlet(double dh, double *phi, int N, double *RHSvec){


    double cc1 = -c/6.0;
    double cc2 = -b/4.0;
    double cc3 = -a/2.0;
    double cc4 =  a/2.0;
    double cc5 =  b/4.0;
    double cc6 =  c/6.0;

    RHSvec[0] = a1*phi[0] + b1*phi[1] + c1*phi[2] + \
                d1*phi[3] + e1*phi[4] + f1*phi[5] + \
		g1*phi[6] + h1*phi[7] + i1*phi[8];

    RHSvec[1] = a2*phi[0] + b2*phi[1] + c2*phi[2] + \
                d2*phi[3] + e2*phi[4] + f2*phi[5] + \
		g2*phi[6] + h2*phi[7];

    RHSvec[2] = a3*phi[0] + b3*phi[1] + c3*phi[2] + \
	        d3*phi[3] + e3*phi[4] + f3*phi[5] + \
		g3*phi[6];

    for(int ip = 3; ip < N-3; ip++){
        RHSvec[ip] = cc1*phi[ip-3] + cc2*phi[ip-2] + cc3*phi[ip-1] + \
                     cc4*phi[ip+1] + cc5*phi[ip+2] + cc6*phi[ip+3];
    }

    RHSvec[N-3] = -g3*phi[N-7] + \
	   	  -f3*phi[N-6] - e3*phi[N-5] - d3*phi[N-4] + \
                  -c3*phi[N-3] - b3*phi[N-2] - a3*phi[N-1];
    RHSvec[N-2] = -h2*phi[N-8] - g2*phi[N-7] + \
	   	  -f2*phi[N-6] - e2*phi[N-5] - d2*phi[N-4] + \
                  -c2*phi[N-3] - b2*phi[N-2] - a2*phi[N-1];
    RHSvec[N-1] = -i1*phi[N-9] - h1*phi[N-8] - g1*phi[N-7] + \
	   	  -f1*phi[N-6] - e1*phi[N-5] - d1*phi[N-4] + \
                  -c1*phi[N-3] - b1*phi[N-2] - a1*phi[N-1];

    for(int ip = 0; ip < N; ip++){
        RHSvec[ip] /= dh;
    }

}

void Penta10::Compact1stPeriodic(double *phi, double *dphidy){

    double RHSvec[N];
    multRHS1stDerivPeriodic(dd, phi, N, RHSvec);

    cyclicPenta(offlower2_1D, offlower_1D, diag_1D, offupper_1D, offupper2_1D, RHSvec, cpvec, dphidy, N);

}

void Penta10::Compact1stPeriodic_TPB(double *phi, double *dphidy){

    double RHSvec[N];
    multRHS1stDerivPeriodic_TPB(dd, phi, N, RHSvec);

    cyclicPenta(offlower2_1D, offlower_1D, diag_1D, offupper_1D, offupper2_1D, RHSvec, cpvec, dphidy, N);

}

void Penta10::Compact1stDirichlet(double *phi, double *dphidy){

    double RHSvec[N];
    double work1[N];
    double work2[N];
    double work3[N];

    multRHS1stDerivDirichlet(dd, phi, N, RHSvec);
    solvePenta(offlower2_1D, offlower_1D, diag_1D, offupper_1D, offupper2_1D, RHSvec, dphidy, work1, work2, work3, N);

}

double Penta10::calcNeumann(double *f){
    return (f[0]*360.0 - f[1]*450.0 + f[2]*400.0 - f[3]*225.0 + f[4]*72.0 - f[5]*10.0)/147.0;   
}
