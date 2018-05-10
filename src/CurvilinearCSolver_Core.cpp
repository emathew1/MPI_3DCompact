#include "CurvilinearCSolver.hpp"

void CurvilinearCSolver::initializeSolverData(){

   
    if(useTiming) ft1 = MPI_Wtime();

    IF_RANK0{
        cout << endl;
        cout << " > Allocating Solver Arrays..." << endl;
        double workSize = 0;
        workSize = 104.0 * (double)N * 8.0;
        cout << " > Need " << workSize/1024.0/1024.0/1024.0 << " Gb of memory required to allocate solver arrays " << endl;
    }

    //3
    c2d->allocY(dU1);
    c2d->allocY(dU2);
    c2d->allocY(dU3);

    //6
    c2d->allocY(dV1);
    c2d->allocY(dV2);
    c2d->allocY(dV3);

    //9
    c2d->allocY(dW1);
    c2d->allocY(dW2);
    c2d->allocY(dW3);

    //12
    c2d->allocY(dT1);
    c2d->allocY(dT2);
    c2d->allocY(dT3);

    //18
    c2d->allocY(Tau11);
    c2d->allocY(Tau12);
    c2d->allocY(Tau13);
    c2d->allocY(Tau22);
    c2d->allocY(Tau23);
    c2d->allocY(Tau33);

    Tau21 = Tau12;
    Tau31 = Tau13;
    Tau32 = Tau23;
   
    //21
    c2d->allocY(cont_1);
    c2d->allocY(cont_2);
    c2d->allocY(cont_3);

    //24
    c2d->allocY(momX_1);
    c2d->allocY(momX_2);
    c2d->allocY(momX_3);

    //27
    c2d->allocY(momY_1);
    c2d->allocY(momY_2);
    c2d->allocY(momY_3);

    //30
    c2d->allocY(momZ_1);
    c2d->allocY(momZ_2);
    c2d->allocY(momZ_3);

    //33
    c2d->allocY(engy_1);
    c2d->allocY(engy_2);
    c2d->allocY(engy_3);

    //37
    c2d->allocY(rho1);
    c2d->allocY(rhok);
    c2d->allocY(rhok2);
    c2d->allocY(rho2);

    //41
    c2d->allocY(rhoU1);
    c2d->allocY(rhoUk);
    c2d->allocY(rhoUk2);
    c2d->allocY(rhoU2);

    //45
    c2d->allocY(rhoV1);
    c2d->allocY(rhoVk);
    c2d->allocY(rhoVk2);
    c2d->allocY(rhoV2);

    //49
    c2d->allocY(rhoW1);
    c2d->allocY(rhoWk);
    c2d->allocY(rhoWk2);
    c2d->allocY(rhoW2);
 
    //53
    c2d->allocY(rhoE1);
    c2d->allocY(rhoEk);
    c2d->allocY(rhoEk2);
    c2d->allocY(rhoE2);

    //58 these will be cleared though...
    c2d->allocY(rho0);
    c2d->allocY(U0);
    c2d->allocY(V0);
    c2d->allocY(W0);
    c2d->allocY(p0);

    //61
    c2d->allocX(U_xp);
    c2d->allocY(U);
    c2d->allocZ(U_zp);

    //64
    c2d->allocX(V_xp);
    c2d->allocY(V);
    c2d->allocZ(V_zp);

    //67
    c2d->allocX(W_xp);
    c2d->allocY(W);
    c2d->allocZ(W_zp);

    //70
    c2d->allocX(T_xp);
    c2d->allocY(T);
    c2d->allocZ(T_zp);

    c2d->allocY(Ucurv);
    c2d->allocY(Vcurv);
    c2d->allocY(Wcurv);

    //72
    c2d->allocY(p);
    c2d->allocY(mu);
    c2d->allocY(sos);

    //84
    c2d->allocX(tempX1);
    c2d->allocX(tempX2);
    c2d->allocX(tempX3);
    c2d->allocX(tempX4);
    c2d->allocX(tempX5);
    c2d->allocX(tempX6);
    c2d->allocX(tempX7);
    c2d->allocX(tempX8);
    c2d->allocX(tempX9);
    c2d->allocX(tempX10);

    //94
    c2d->allocY(tempY1);
    c2d->allocY(tempY2);
    c2d->allocY(tempY3);
    c2d->allocY(tempY4);
    c2d->allocY(tempY5);
    c2d->allocY(tempY6);
    c2d->allocY(tempY7);
    c2d->allocY(tempY8);
    c2d->allocY(tempY9);
    c2d->allocY(tempY10);
    c2d->allocY(tempY11);
    c2d->allocY(tempY12);

    //104
    c2d->allocZ(tempZ1);
    c2d->allocZ(tempZ2);
    c2d->allocZ(tempZ3);
    c2d->allocZ(tempZ4);
    c2d->allocZ(tempZ5);
    c2d->allocZ(tempZ6);
    c2d->allocZ(tempZ7);
    c2d->allocZ(tempZ8);
    c2d->allocZ(tempZ9);
    c2d->allocZ(tempZ10);

    
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > initSolDat Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
}


void CurvilinearCSolver::calcDtFromCFL(){

    if(useTiming) ft1 = MPI_Wtime();    

    //Calculate the wave speed over the local spacings...
    double *UChar_dx;
    c2d->allocY(UChar_dx);

    FOR_XYZ_YPEN{
	UChar_dx[ip] = (fabs(U[ip]) + sos[ip])/dom->dx + (fabs(V[ip])+sos[ip])/dom->dy + (fabs(W[ip]) + sos[ip])/dom->dz;
    }

    //Get the largest value in the domain
    double lmax_UChar_dx = -100000.0;
    FOR_XYZ_YPEN{
	if(UChar_dx[ip] > lmax_UChar_dx){
	    lmax_UChar_dx = UChar_dx[ip];
	}
    }
    
    double max_UChar_dx;
    MPI_Allreduce(&lmax_UChar_dx, &max_UChar_dx, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    //done with UChar_dx
    c2d->deallocXYZ(UChar_dx);

    if(ts->CONST_CFL){
	ts->dt = ts->CFL/max_UChar_dx;
    }else if(ts->CONST_DT){
	ts->CFL = ts->dt*max_UChar_dx;
    }
  
    if(timeStep == 0){
	timeStep++;
	time = 0.0;
	IF_RANK0 cout << endl;
    }else{
	timeStep++;
	time += ts->dt;
    }

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > calcDtCFL  Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }



}

void CurvilinearCSolver::preStepDerivatives(){

    if(useTiming) ft1 = MPI_Wtime();

    double *rhoP;
    double *rhoUP;
    double *rhoVP;
    double *rhoWP;
    double *rhoEP;

    if(rkStep == 1){
	rhoP  = rho1;
	rhoUP = rhoU1;
	rhoVP = rhoV1;
	rhoWP = rhoW1;
	rhoEP = rhoE1; 
    }else{
	rhoP  = rhok;
	rhoUP = rhoUk;
	rhoVP = rhoVk;
	rhoWP = rhoWk;
	rhoEP = rhoEk; 
    }

    double ftt1 = 0, ftt2 = 0;

    if(useTiming) ftt1 = MPI_Wtime();

    ///////////////////
    // Y-DERIVATIVES //
    ///////////////////


    //First we'll do all of the Y-Direction derivatives to calc tau
    derivY->calc1stDerivField(U, Uy); //dU/dx
    derivY->calc1stDerivField(V, Vy); //dV/dx
    derivY->calc1stDerivField(W, Wy); //dW/dx
    derivY->calc1stDerivField(T, Ty); //dT/dx


    if(useTiming){
	ftt2 = MPI_Wtime();
	IF_RANK0 cout << " > yderivtrans Timing: " << setw(6)  << (int)((ftt2-ftt1)*1000) << "ms" << endl;
	ftt1 = MPI_Wtime();
    }

    ///////////////////
    // X-DERIVATIVES //
    ///////////////////

    double *dU11, *Vx1, *Wx1, *Tx1;
    double *U1, *V1, *W1, *T1;

    //Point to the needed X memory
    U1 = tempX1; Ux1 = tempX5; 
    V1 = tempX2; Vx1 = tempX6;
    W1 = tempX3; Wx1 = tempX7;
    T1 = tempX4; Tx1 = tempX8;

    c2d->transposeY2X_MajorIndex(U, U1);
    c2d->transposeY2X_MajorIndex(V, V1);
    c2d->transposeY2X_MajorIndex(W, W1);
    c2d->transposeY2X_MajorIndex(T, T1);

    derivX->calc1stDerivField(U1, Ux1);
    derivX->calc1stDerivField(V1, Vx1);
    derivX->calc1stDerivField(W1, Wx1);
    derivX->calc1stDerivField(T1, Tx1);

    c2d->transposeX2Y_MajorIndex(Ux1, Ux);
    c2d->transposeX2Y_MajorIndex(Vx1, Vx);
    c2d->transposeX2Y_MajorIndex(Wx1, Wx);
    c2d->transposeX2Y_MajorIndex(Tx1, Tx);

    if(useTiming){
	ftt2 = MPI_Wtime();
	IF_RANK0 cout << " > xderivtrans Timing: " << setw(6)  << (int)((ftt2-ftt1)*1000) << "ms" << endl;
	ftt1 = MPI_Wtime();
    }


    ///////////////////
    // Z-DERIVATIVES //
    ///////////////////

    double *Uz3, *Vz3, *Wz3, *Tz3;
    double *U3,  *V3,  *W3,  *T3;

    //Point to the needed Z memory
    Uz3 = tempZ1;
    Vz3 = tempZ2;
    Wz3 = tempZ3;
    Tz3 = tempZ4;

    //Point to the Z memory
    U3 = tempZ5;
    V3 = tempZ6;
    W3 = tempZ7;
    T3 = tempZ8; 

    c2d->transposeY2Z_MajorIndex(U, U3);
    c2d->transposeY2Z_MajorIndex(V, V3);
    c2d->transposeY2Z_MajorIndex(W, W3);
    c2d->transposeY2Z_MajorIndex(T, T3);

    derivZ->calc1stDerivField(U3, Uz3);
    derivZ->calc1stDerivField(V3, Vz3);
    derivZ->calc1stDerivField(W3, Wz3);
    derivZ->calc1stDerivField(T3, Tz3);

    c2d->transposeZ2Y_MajorIndex(Uz3, Uz);
    c2d->transposeZ2Y_MajorIndex(Vz3, Vz);
    c2d->transposeZ2Y_MajorIndex(Wz3, Wz);
    c2d->transposeZ2Y_MajorIndex(Tz3, Tz);

    if(useTiming){
	ftt2 = MPI_Wtime();
	IF_RANK0 cout << " > zderivtrans Timing: " << setw(6)  << (int)((ftt2-ftt1)*1000) << "ms" << endl;
	ftt1 = MPI_Wtime();
    }


    /////////////////////////////////////////
    // CALC TAU COMPONENTS & EQN COMPONENTS//
    /////////////////////////////////////////

    double *preMomX_X, *preMomX_Y, *preMomX_Z; 
    double *preMomY_X, *preMomY_Y, *preMomY_Z; 
    double *preMomZ_X, *preMomZ_Y, *preMomZ_Z; 
    double *preEngy_X, *preEngy_Y, *preEngy_Z; 

    preMomX_X = tempY1;  preMomX_Y = tempY2;  preMomX_Z = tempY3;
    preMomY_X = tempY4;  preMomY_Y = tempY5;  preMomY_Z = tempY6;
    preMomZ_X = tempY7;  preMomZ_Y = tempY8;  preMomZ_Z = tempY9;
    preEngy_X = tempY10; preEngy_Y = tempY11; preEngy_Z = tempY12;

    //Now recalculate properties in the new space
    FOR_XYZ_YPEN{
	Tauxx[ip] = mu[ip]*((4.0/3.0)*Ux[ip] - (2.0/3.0)*(Vy[ip] + Wz[ip]));
	Tauyy[ip] = mu[ip]*((4.0/3.0)*Vy[ip] - (2.0/3.0)*(Ux[ip] + Wz[ip]));
	Tauzz[ip] = mu[ip]*((4.0/3.0)*Wz[ip] - (2.0/3.0)*(Ux[ip] + Vy[ip]));
	Tauxy[ip] = mu[ip]*(Uy[ip] + Vx[ip]);
	Tauxz[ip] = mu[ip]*(Uz[ip] + Wx[ip]);
	Tauyz[ip] = mu[ip]*(Vz[ip] + Wy[ip]);

	preMomX_X[ip] = rhoUP[ip]*U[ip] + p[ip] - Tauxx[ip];
	preMomX_Y[ip] = rhoUP[ip]*V[ip] - Tauxy[ip];
	preMomX_Z[ip] = rhoUP[ip]*W[ip] - Tauxz[ip];

	preMomY_X[ip] = rhoVP[ip]*U[ip] - Tauyx[ip];
	preMomY_Y[ip] = rhoVP[ip]*V[ip] + p[ip] - Tauyy[ip];
	preMomY_Z[ip] = rhoVP[ip]*W[ip] - Tauyz[ip];

	preMomZ_X[ip] = rhoWP[ip]*U[ip] - Tauzx[ip];
	preMomZ_Y[ip] = rhoWP[ip]*V[ip] - Tauzy[ip];
	preMomZ_Z[ip] = rhoWP[ip]*W[ip] + p[ip] - Tauzz[ip];

	preEngy_X[ip] = rhoEP[ip]*U[ip] + U[ip]*p[ip] - (ig->cp/ig->Pr)*mu[ip]*Tx[ip] - U[ip]*Tauxx[ip] - V[ip]*Tauyx[ip] - W[ip]*Tauzx[ip];
	preEngy_Y[ip] = rhoEP[ip]*V[ip] + V[ip]*p[ip] - (ig->cp/ig->Pr)*mu[ip]*Ty[ip] - U[ip]*Tauxy[ip] - V[ip]*Tauyy[ip] - W[ip]*Tauzy[ip];
	preEngy_Z[ip] = rhoEP[ip]*W[ip] + W[ip]*p[ip] - (ig->cp/ig->Pr)*mu[ip]*Tz[ip] - U[ip]*Tauxz[ip] - V[ip]*Tauyz[ip] - W[ip]*Tauzz[ip];

    }

    if(useTiming){
	ftt2 = MPI_Wtime();
	IF_RANK0 cout << " > Tau&Components Timing: " << setw(6)  << (int)((ftt2-ftt1)*1000) << "ms" << endl;
    }

    ///////////////////
    // Y-DERIVATIVES //
    ///////////////////
   
    derivY->calc1stDerivField(rhoVP,     cont_Y);
    derivY->calc1stDerivField(preMomX_Y, momX_Y);
    derivY->calc1stDerivField(preMomY_Y, momY_Y);
    derivY->calc1stDerivField(preMomZ_Y, momZ_Y);
    derivY->calc1stDerivField(preEngy_Y, engy_Y);

    if(useTiming){
	ftt2 = MPI_Wtime();
	IF_RANK0 cout << " > yderivtrans Timing: " << setw(6)  << (int)((ftt2-ftt1)*1000) << "ms" << endl;
    }

    ///////////////////
    // X-DERIVATIVES //
    ///////////////////

    double *rhoUP1, *preMomX_X1, *preMomY_X1, *preMomZ_X1, *preEngy_X1;
    rhoUP1     = tempX1; 
    preMomX_X1 = tempX2; 
    preMomY_X1 = tempX3; 
    preMomZ_X1 = tempX4; 
    preEngy_X1 = tempX5; 

    double *cont_X1, *momX_X1, *momY_X1, *momZ_X1, *engy_X1;
    cont_X1 = tempX6; 
    momX_X1 = tempX7; 
    momY_X1 = tempX8; 
    momZ_X1 = tempX9; 
    engy_X1 = tempX10; 

    c2d->transposeY2X_MajorIndex(rhoUP,     rhoUP1);
    c2d->transposeY2X_MajorIndex(preMomX_X, preMomX_X1);
    c2d->transposeY2X_MajorIndex(preMomY_X, preMomY_X1);
    c2d->transposeY2X_MajorIndex(preMomZ_X, preMomZ_X1);
    c2d->transposeY2X_MajorIndex(preEngy_X, preEngy_X1);

    derivX->calc1stDerivField(rhoUP1,     cont_X1); 
    derivX->calc1stDerivField(preMomX_X1, momX_X1); 
    derivX->calc1stDerivField(preMomY_X1, momY_X1); 
    derivX->calc1stDerivField(preMomZ_X1, momZ_X1); 
    derivX->calc1stDerivField(preEngy_X1, engy_X1); 

    c2d->transposeX2Y_MajorIndex(cont_X1, cont_X);
    c2d->transposeX2Y_MajorIndex(momX_X1, momX_X);
    c2d->transposeX2Y_MajorIndex(momY_X1, momY_X);
    c2d->transposeX2Y_MajorIndex(momZ_X1, momZ_X);
    c2d->transposeX2Y_MajorIndex(engy_X1, engy_X);

    if(useTiming){
	ftt2 = MPI_Wtime();
	IF_RANK0 cout << " > xderivtrans Timing: " << setw(6)  << (int)((ftt2-ftt1)*1000) << "ms" << endl;
    }

    ///////////////////
    // Z-DERIVATIVES //
    ///////////////////

    double *rhoWP3, *preMomX_Z3, *preMomY_Z3, *preMomZ_Z3, *preEngy_Z3;
    double *cont_Z3, *momX_Z3, *momY_Z3, *momZ_Z3, *engy_Z3;

    rhoWP3     = tempZ1; cont_Z3 = tempZ6;
    preMomX_Z3 = tempZ2; momX_Z3 = tempZ7;
    preMomY_Z3 = tempZ3; momY_Z3 = tempZ8;
    preMomZ_Z3 = tempZ4; momZ_Z3 = tempZ9;
    preEngy_Z3 = tempZ5; engy_Z3 = tempZ10;

    c2d->transposeY2Z_MajorIndex(rhoWP,     rhoWP3);
    c2d->transposeY2Z_MajorIndex(preMomX_Z, preMomX_Z3);
    c2d->transposeY2Z_MajorIndex(preMomY_Z, preMomY_Z3);
    c2d->transposeY2Z_MajorIndex(preMomZ_Z, preMomZ_Z3);
    c2d->transposeY2Z_MajorIndex(preEngy_Z, preEngy_Z3);

    derivZ->calc1stDerivField(rhoWP3,     cont_Z3);
    derivZ->calc1stDerivField(preMomX_Z3, momX_Z3);
    derivZ->calc1stDerivField(preMomY_Z3, momY_Z3);
    derivZ->calc1stDerivField(preMomZ_Z3, momZ_Z3);
    derivZ->calc1stDerivField(preEngy_Z3, engy_Z3);

    c2d->transposeZ2Y_MajorIndex(cont_Z3, cont_Z);
    c2d->transposeZ2Y_MajorIndex(momX_Z3, momX_Z);
    c2d->transposeZ2Y_MajorIndex(momY_Z3, momY_Z);
    c2d->transposeZ2Y_MajorIndex(momZ_Z3, momZ_Z);
    c2d->transposeZ2Y_MajorIndex(engy_Z3, engy_Z);

    if(useTiming){
	ftt2 = MPI_Wtime();
	IF_RANK0 cout << " > zderivtrans Timing: " << setw(6)  << (int)((ftt2-ftt1)*1000) << "ms" << endl;
    }


   

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > preStepDer Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }


}


void CurvilinearCSolver::solveContinuity(){

    if(useTiming) ft1 = MPI_Wtime();

    double *rhoP;
    if(rkStep == 1){
        rhoP = rho1;
    }else{
        rhoP = rhok;
    }
	
    double spgSource; 

    FOR_XYZ_YPEN{

	if(spongeFlag)
	    spgSource = calcSpongeSource(rhoP[ip], spg->spongeRhoAvg[ip], spg->sigma[ip]);
	else
	    spgSource = 0.0;
		
	rhok2[ip]  = ts->dt*(-cont_X[ip] - cont_Y[ip] - cont_Z[ip] + spgSource);
    }

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > solveCont  Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }

}

void CurvilinearCSolver::solveXMomentum(){

    if(useTiming) ft1 = MPI_Wtime();

    double *rhoUP;
    if(rkStep == 1){
        rhoUP = rhoU1;
    }else{
        rhoUP = rhoUk;
    }

    double spgSource;
    FOR_XYZ_YPEN{

	if(spongeFlag)
	    spgSource = calcSpongeSource(rhoUP[ip], spg->spongeRhoUAvg[ip], spg->sigma[ip]);
	else
	    spgSource = 0.0;

	rhoUk2[ip] += -momX_X[ip] - momX_Y[ip] -momX_Z[ip] + spgSource;
	rhoUk2[ip] *= ts->dt;
    }

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > solveXMom  Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }



}

void CurvilinearCSolver::solveYMomentum(){

    if(useTiming) ft1 = MPI_Wtime();

    double *rhoVP;
    if(rkStep == 1){
        rhoVP = rhoV1;
    }else{
        rhoVP = rhoVk;
    }

    double spgSource;

    FOR_XYZ_YPEN{ 

        if(spongeFlag)
            spgSource = calcSpongeSource(rhoVP[ip], spg->spongeRhoVAvg[ip], spg->sigma[ip]);
	else
	    spgSource = 0.0;

	rhoVk2[ip] += -momY_X[ip] -momY_Y[ip] -momY_Z[ip] + spgSource;
        rhoVk2[ip] *= ts->dt;
   }

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > solveYMom  Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }



}

void CurvilinearCSolver::solveZMomentum(){

    if(useTiming) ft1 = MPI_Wtime();

    double *rhoWP;
    if(rkStep == 1){
        rhoWP = rhoW1;
    }else{
        rhoWP = rhoWk;
    }

    double spgSource;
    FOR_XYZ_YPEN{

        if(spongeFlag)
            spgSource = calcSpongeSource(rhoWP[ip], spg->spongeRhoWAvg[ip], spg->sigma[ip]);
        else
	    spgSource = 0.0;

	rhoWk2[ip] += -momZ_X[ip] -momZ_Y[ip] -momZ_Z[ip] + spgSource;
        rhoWk2[ip] *= ts->dt;
    }

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > solveZMom  Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }



}


void CurvilinearCSolver::solveEnergy(){

    if(useTiming) ft1 = MPI_Wtime();

    double *rhoEP;
    if(rkStep == 1){
        rhoEP = rhoE1;
    }else{
        rhoEP = rhoEk;
    }

    double spgSource;
    FOR_XYZ_YPEN{

        if(spongeFlag)
            spgSource = calcSpongeSource(rhoEP[ip], spg->spongeRhoEAvg[ip], spg->sigma[ip]);
        else
            spgSource = 0.0;

	rhoEk2[ip] = -engy_X[ip] - engy_Y[ip] - engy_Z[ip] + spgSource;
	rhoEk2[ip] *= ts->dt;
    }

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > solveEngy  Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }


}

void CurvilinearCSolver::filterConservedData(){

    if(useTiming) ft1 = MPI_Wtime();

    //Need to do round robin of filtering directions...
    if(timeStep%ts->filterStep == 0){

        //Advance the filtering time step       
        filterTimeStep++;

        //Going to try and be cute to minimize dmemory allocation
        if(filterTimeStep%3 == 1){

            //Here we'll do Y->Z->X     

                    filtY->filterField(rho2,  tempY1);
		    filtY->filterField(rhoU2, tempY2);
		    filtY->filterField(rhoV2, tempY3);
		    filtY->filterField(rhoW2, tempY4);
		    filtY->filterField(rhoE2, tempY5);

		    c2d->transposeY2Z_MajorIndex(tempY1, tempZ1);
		    c2d->transposeY2Z_MajorIndex(tempY2, tempZ2);
		    c2d->transposeY2Z_MajorIndex(tempY3, tempZ3);
		    c2d->transposeY2Z_MajorIndex(tempY4, tempZ4);
		    c2d->transposeY2Z_MajorIndex(tempY5, tempZ5);

		    filtZ->filterField(tempZ1, tempZ6);
		    filtZ->filterField(tempZ2, tempZ7);
		    filtZ->filterField(tempZ3, tempZ8);
		    filtZ->filterField(tempZ4, tempZ9);
		    filtZ->filterField(tempZ5, tempZ10);

		    c2d->transposeZ2Y_MajorIndex(tempZ6,  tempY1);
		    c2d->transposeZ2Y_MajorIndex(tempZ7,  tempY2);
		    c2d->transposeZ2Y_MajorIndex(tempZ8,  tempY3);
		    c2d->transposeZ2Y_MajorIndex(tempZ9,  tempY4);
		    c2d->transposeZ2Y_MajorIndex(tempZ10, tempY5);

		    c2d->transposeY2X_MajorIndex(tempY1, tempX1);
		    c2d->transposeY2X_MajorIndex(tempY2, tempX2);
		    c2d->transposeY2X_MajorIndex(tempY3, tempX3);
		    c2d->transposeY2X_MajorIndex(tempY4, tempX4);
		    c2d->transposeY2X_MajorIndex(tempY5, tempX5);
		  
	            filtX->filterField(tempX1, tempX6);	
	            filtX->filterField(tempX2, tempX7);	
	            filtX->filterField(tempX3, tempX8);	
	            filtX->filterField(tempX4, tempX9);	
	            filtX->filterField(tempX5, tempX10);	

		    c2d->transposeX2Y_MajorIndex(tempX6,  rho1);
		    c2d->transposeX2Y_MajorIndex(tempX7,  rhoU1);
		    c2d->transposeX2Y_MajorIndex(tempX8,  rhoV1);
		    c2d->transposeX2Y_MajorIndex(tempX9,  rhoW1);
		    c2d->transposeX2Y_MajorIndex(tempX10, rhoE1);


        }else if(filterTimeStep%3 == 2){

            //Here we'll do Z->X->Y     

	    c2d->transposeY2Z_MajorIndex(rho2,  tempZ1);
	    c2d->transposeY2Z_MajorIndex(rhoU2, tempZ2);
	    c2d->transposeY2Z_MajorIndex(rhoV2, tempZ3);
	    c2d->transposeY2Z_MajorIndex(rhoW2, tempZ4);
	    c2d->transposeY2Z_MajorIndex(rhoE2, tempZ5);

	    filtZ->filterField(tempZ1, tempZ6);
	    filtZ->filterField(tempZ2, tempZ7);
	    filtZ->filterField(tempZ3, tempZ8);
	    filtZ->filterField(tempZ4, tempZ9);
	    filtZ->filterField(tempZ5, tempZ10);

	    c2d->transposeZ2Y_MajorIndex(tempZ6,  tempY1);
	    c2d->transposeZ2Y_MajorIndex(tempZ7,  tempY2);
	    c2d->transposeZ2Y_MajorIndex(tempZ8,  tempY3);
	    c2d->transposeZ2Y_MajorIndex(tempZ9,  tempY4);
	    c2d->transposeZ2Y_MajorIndex(tempZ10, tempY5);

	    c2d->transposeY2X_MajorIndex(tempY1,  tempX1);
	    c2d->transposeY2X_MajorIndex(tempY2,  tempX2);
	    c2d->transposeY2X_MajorIndex(tempY3,  tempX3);
	    c2d->transposeY2X_MajorIndex(tempY4,  tempX4);
	    c2d->transposeY2X_MajorIndex(tempY5,  tempX5);

	    filtX->filterField(tempX1, tempX6);
	    filtX->filterField(tempX2, tempX7);
	    filtX->filterField(tempX3, tempX8);
	    filtX->filterField(tempX4, tempX9);
	    filtX->filterField(tempX5, tempX10);

	    c2d->transposeX2Y_MajorIndex(tempX6,  tempY1);
	    c2d->transposeX2Y_MajorIndex(tempX7,  tempY2);
	    c2d->transposeX2Y_MajorIndex(tempX8,  tempY3);
	    c2d->transposeX2Y_MajorIndex(tempX9,  tempY4);
	    c2d->transposeX2Y_MajorIndex(tempX10, tempY5);

	    filtY->filterField(tempY1,  rho1);
	    filtY->filterField(tempY2, rhoU1);
	    filtY->filterField(tempY3, rhoV1);
	    filtY->filterField(tempY4, rhoW1);
	    filtY->filterField(tempY5, rhoE1);

        }else{

            //Here we'll do X->Y->Z     

	    c2d->transposeY2X_MajorIndex(rho2,  tempX1);
	    c2d->transposeY2X_MajorIndex(rhoU2, tempX2);
	    c2d->transposeY2X_MajorIndex(rhoV2, tempX3);
	    c2d->transposeY2X_MajorIndex(rhoW2, tempX4);
	    c2d->transposeY2X_MajorIndex(rhoE2, tempX5);

	    filtX->filterField(tempX1, tempX6);
	    filtX->filterField(tempX2, tempX7);
	    filtX->filterField(tempX3, tempX8);
	    filtX->filterField(tempX4, tempX9);
	    filtX->filterField(tempX5, tempX10);

	    c2d->transposeX2Y_MajorIndex(tempX6,  tempY1);
	    c2d->transposeX2Y_MajorIndex(tempX7,  tempY2);
	    c2d->transposeX2Y_MajorIndex(tempX8,  tempY3);
	    c2d->transposeX2Y_MajorIndex(tempX9,  tempY4);
	    c2d->transposeX2Y_MajorIndex(tempX10, tempY5);

	    filtY->filterField(tempY1, tempY6);
	    filtY->filterField(tempY2, tempY7);
	    filtY->filterField(tempY3, tempY8);
	    filtY->filterField(tempY4, tempY9);
	    filtY->filterField(tempY5, tempY10);

	    c2d->transposeY2Z_MajorIndex(tempY6,  tempZ1);
	    c2d->transposeY2Z_MajorIndex(tempY7,  tempZ2);
	    c2d->transposeY2Z_MajorIndex(tempY8,  tempZ3);
	    c2d->transposeY2Z_MajorIndex(tempY9,  tempZ4);
	    c2d->transposeY2Z_MajorIndex(tempY10, tempZ5);

	    filtZ->filterField(tempZ1, tempZ6);
	    filtZ->filterField(tempZ2, tempZ7);
	    filtZ->filterField(tempZ3, tempZ8);
	    filtZ->filterField(tempZ4, tempZ9);
	    filtZ->filterField(tempZ5, tempZ10);

	    c2d->transposeZ2Y_MajorIndex(tempZ6,  rho1);
	    c2d->transposeZ2Y_MajorIndex(tempZ7,  rhoU1);
	    c2d->transposeZ2Y_MajorIndex(tempZ8,  rhoV1);
	    c2d->transposeZ2Y_MajorIndex(tempZ9,  rhoW1);
	    c2d->transposeZ2Y_MajorIndex(tempZ10, rhoE1);

        }
 
    //If not filtering, need to copy the solution over to the *1 variables
    }else{
	    memcpy(rho1,  rho2,  sizeof(double)*pySize[0]*pySize[1]*pySize[2]);
	    memcpy(rhoU1, rhoU2, sizeof(double)*pySize[0]*pySize[1]*pySize[2]);
	    memcpy(rhoV1, rhoV2, sizeof(double)*pySize[0]*pySize[1]*pySize[2]);
	    memcpy(rhoW1, rhoW2, sizeof(double)*pySize[0]*pySize[1]*pySize[2]);
	    memcpy(rhoE1, rhoE2, sizeof(double)*pySize[0]*pySize[1]*pySize[2]);

    }

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > filterCons Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }



};


void CurvilinearCSolver::updateNonConservedData(){

    if(useTiming) ft1 = MPI_Wtime();

    if(!rkLast){

	FOR_XYZ_YPEN{
	    U[ip]   = ig->solveU(rhok[ip], rhoUk[ip]);
	    V[ip]   = ig->solveU(rhok[ip], rhoVk[ip]);
	    W[ip]   = ig->solveU(rhok[ip], rhoWk[ip]);
	    p[ip]   = ig->solvep(rhok[ip], rhoEk[ip], U[ip], V[ip], W[ip]);
	    T[ip]   = ig->solveT(rhok[ip], p[ip]);
	    mu[ip]  = ig->solveMu(T[ip]);
	    sos[ip] = ig->solveSOS(rhok[ip], p[ip]);
	}

    }else if(rkLast){

	FOR_XYZ_YPEN{
	    U[ip]   = ig->solveU(rho1[ip], rhoU1[ip]);
	    V[ip]   = ig->solveU(rho1[ip], rhoV1[ip]);
	    W[ip]   = ig->solveU(rho1[ip], rhoW1[ip]);
	    p[ip]   = ig->solvep(rho1[ip], rhoE1[ip], U[ip], V[ip], W[ip]);
	    T[ip]   = ig->solveT(rho1[ip], p[ip]);
	    mu[ip]  = ig->solveMu(T[ip]);
	    sos[ip] = ig->solveSOS(rho1[ip], p[ip]);
	}
    }

    
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > updNonCons Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }



}

void CurvilinearCSolver::checkSolution(){


    if(useTiming) ft1 = MPI_Wtime();

    if(timeStep%ts->checkStep == 0){

	t2 = MPI_Wtime();

	IF_RANK0{
	    cout << endl;
            cout << "-------------------------------------------------" << endl;
            cout << " Step = "<< timeStep << ", time = " << time << ", dt = " << ts->dt << endl;
            cout << "-------------------------------------------------" << endl;
            cout << "  Time since last timestep = " << t2 - t1  << endl;
	}
	

	int Nx = pySize[0];
	int Ny = pySize[1];
	int Nz = pySize[2];

        getRange(rho1, "RHO", Nx, Ny, Nz, mpiRank);
        getRange(U, "U", Nx, Ny, Nz, mpiRank);
        getRange(V, "V", Nx, Ny, Nz, mpiRank);
        getRange(W, "W", Nx, Ny, Nz, mpiRank);
        getRange(p, "P", Nx, Ny, Nz, mpiRank);
        getRange(T, "T", Nx, Ny, Nz, mpiRank);
        getRange(mu, "mu", Nx, Ny, Nz, mpiRank);
        getRange(rhoE1, "RHOE", Nx, Ny, Nz, mpiRank);
        getRange(sos, "SOS", Nx, Ny, Nz, mpiRank);
        IF_RANK0 cout << endl;

	t1 = MPI_Wtime();
    }

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > checkSoln  Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }	 


};


void CurvilinearCSolver::dumpSolution(){

    if(useTiming) ft1 = MPI_Wtime();

	double time1 = MPI_Wtime();

	IF_RANK0{
            cout << endl;
            cout << " > ===============" << endl;
            cout << " >  DUMPING FIELD " << endl;
            cout << " > ===============" << endl;
	}

        ofstream outfile;
        outfile.precision(17);
        string outputFileName;
	outputFileName = "SolutionDump.";
	
	ostringstream timeStepString;
        timeStepString << timeStep;

	outputFileName.append(timeStepString.str());

	MPI_File fh;
	MPI_Offset disp, filesize;

	MPI_File_open(MPI_COMM_WORLD, outputFileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	filesize = 0;
	MPI_File_set_size(fh, filesize);

	disp = 0;
	c2d->writeVar(fh, disp, baseDirection, rho1);
	c2d->writeVar(fh, disp, baseDirection, rhoU1);
	c2d->writeVar(fh, disp, baseDirection, rhoV1);
	c2d->writeVar(fh, disp, baseDirection, rhoW1);
	c2d->writeVar(fh, disp, baseDirection, rhoE1);

	MPI_File_close(&fh);

	double time2 = MPI_Wtime();

	IF_RANK0{
	    cout << endl;
	    cout << " > File dump took " << time2-time1 << " seconds." << endl;
	}

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > dumpSoln Timing: " << setw(6)  << endl;
	
    }	 



}

void CurvilinearCSolver::writeImages(){

    if(useTiming) ft1 = MPI_Wtime();

	IF_RANK0 cout << " > Dumping images..." << endl;

	ostringstream timeStepStringT;
	timeStepStringT << timeStep;
	int zeroPad = 6;
	string timeStepString = string(zeroPad - (timeStepStringT.str()).length(), '0') + timeStepStringT.str();

	writePlaneImageForVariable(p, "pXZ", timeStepString, 1, 0.5);
	writePlaneImageForVariable(U, "uXZ", timeStepString, 1, 0.5);
	writePlaneImageForVariable(p, "pYZ", timeStepString, 0, 0.5);

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > writeImg Timing: " << setw(6)  << endl;
    }	 



}

void CurvilinearCSolver::writePlaneImageForVariable(double *var, string varName, string timeStepString, int plane, double fraction){

    //plane == 0, YZ (X-normal Plane)
    //plane == 1, XZ (Y-normal Plane)
    //plane == 2, XY (Z-normal Plane)

    int N1 = 0, N2 = 0;
    PngWriter *png = NULL;

    if(plane == 0){
	N1 = Ny;
	N2 = Nz;
	
	IF_RANK0 png = pngYZ;	

    }else if(plane == 1){
	N1 = Nz;
	N2 = Nx;

	IF_RANK0 png = pngXZ;	

    }else if(plane == 2){
	N1 = Nx;
	N2 = Ny;

	IF_RANK0 png = pngXY;	
    }
    
    double *ff = new double[N1*N2];
    for(int ip = 0; ip < N1*N2; ip++) ff[ip] = -1000000.0;


    FOR_Z_YPEN{
	FOR_Y_YPEN{
	    FOR_X_YPEN{
		int gi = GETGLOBALXIND_YPEN;
		int gj = GETGLOBALYIND_YPEN;
		int gk = GETGLOBALZIND_YPEN;
	
		int ip = GETMAJIND_YPEN;

		if(plane == 0){
		    if(gi == (int)(Nx*fraction)){
			ff[gk*Ny + gj] = var[ip];
		    }
		}else if(plane == 1){
		    if(gj == (int)(Ny*fraction)){
			ff[gi*Nz + gk] = var[ip];
		    }
		}else if(plane == 2){
		    if(gk == (int)(Nz*fraction)){
			ff[gj*Nx + gi] = var[ip];
		    }
		}else{
		   IF_RANK0 cout << "Unrecognized plane value, must be 1, 2, or 3! May crash? " << endl;
		}

	    }
	}
    }

    double *gf;
    IF_RANK0{
	gf = new double[N1*N2];
    }else{
  	gf = NULL;
    }

    MPI_Reduce(ff, gf, N1*N2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    IF_RANK0{

	//Make sure all of the values were filled
	for(int ip = 0; ip < N1*N2; ip++){

	    if(gf[ip] == -1000000.0){
		cout << "Error: image matrix not completely filled!" << endl;
	    }
	}

	//get the max/min value in the plane...
	double dataMin =  100000000.0;
	double dataMax = -100000000.0;
	for(int ip = 0; ip < N1*N2; ip++){
	    dataMin = fmin(dataMin, gf[ip]);
	    dataMax = fmax(dataMax, gf[ip]);
	}

	//Scale pixel value to local min and max
	int *g = new int[N1*N2];
	for(int ip = 0; ip < N1*N2; ip++){
	   double gtemp = (gf[ip] - dataMin)/(dataMax - dataMin); 
	   g[ip] = (int)(gtemp*255.0);
	}

	for(int jp = 0; jp < N2; jp++){
	    for(int ip = 0; ip < N1; ip++){
		int ii = jp*N1 + ip;
		//Grayscale image...
		png->set(ip, jp, g[ii], g[ii], g[ii]);
	    }
	}

	string imageName = varName;
	imageName.append(".");
	imageName.append(timeStepString);
	imageName.append(".png");
	png->write(imageName.c_str());

	delete[] g;
	delete[] gf;
    }

    delete[] ff;

}

void CurvilinearCSolver::checkEnd(){


    if(time >= ts->maxTime){

	IF_RANK0{
	    cout << "=================" << endl;
	    cout << " HIT END OF TIME " << endl;
	    cout << "=================" << endl;
	}

	endFlag = true;
    }

    if(timeStep >= ts->maxTimeStep){
	
	IF_RANK0{
	    cout << "=================" << endl;
	    cout << " HIT END OF TIME " << endl;
	    cout << "=================" << endl;
	}

	endFlag = true;

    } 

    if(endFlag){
	dumpSolution();
    }

}

void CurvilinearCSolver::reportAll(){

   int Nx = pySize[0];
   int Ny = pySize[1];
   int Nz = pySize[2];

   IF_RANK0   cout << "REPORT ALL" << endl;

   getRange(Ux, "Ux", Nx, Ny, Nz, mpiRank);
   getRange(Uy, "Uy", Nx, Ny, Nz, mpiRank);
   getRange(Uz, "Uz", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(Vx, "Vx", Nx, Ny, Nz, mpiRank);
   getRange(Vy, "Vy", Nx, Ny, Nz, mpiRank);
   getRange(Vz, "Vz", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(Wx, "Wx", Nx, Ny, Nz, mpiRank);
   getRange(Wy, "Wy", Nx, Ny, Nz, mpiRank);
   getRange(Wz, "Wz", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(Tx, "Tx", Nx, Ny, Nz, mpiRank);
   getRange(Ty, "Ty", Nx, Ny, Nz, mpiRank);
   getRange(Tz, "Tz", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(Tauxx, "Tauxx", Nx, Ny, Nz, mpiRank);
   getRange(Tauyy, "Tauyy", Nx, Ny, Nz, mpiRank);
   getRange(Tauzz, "Tauzz", Nx, Ny, Nz, mpiRank);
   getRange(Tauxy, "Tauxy", Nx, Ny, Nz, mpiRank);
   getRange(Tauyz, "Tauyz", Nx, Ny, Nz, mpiRank);
   getRange(Tauxz, "Tauxz", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(cont_X, "cont_X", Nx, Ny, Nz, mpiRank);
   getRange(cont_Y, "cont_Y", Nx, Ny, Nz, mpiRank);
   getRange(cont_Z, "cont_Z", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(momX_X, "momX_X", Nx, Ny, Nz, mpiRank);
   getRange(momX_Y, "momX_Y", Nx, Ny, Nz, mpiRank);
   getRange(momX_Z, "momX_Z", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(momY_X, "momY_X", Nx, Ny, Nz, mpiRank);
   getRange(momY_Y, "momY_Y", Nx, Ny, Nz, mpiRank);
   getRange(momY_Z, "momY_Z", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(momZ_X, "momZ_X", Nx, Ny, Nz, mpiRank);
   getRange(momZ_Y, "momZ_Y", Nx, Ny, Nz, mpiRank);
   getRange(momZ_Z, "momZ_Z", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(engy_X, "engy_X", Nx, Ny, Nz, mpiRank);
   getRange(engy_Y, "engy_Y", Nx, Ny, Nz, mpiRank);
   getRange(engy_Z, "engy_Z", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(rho1, "rho1", Nx, Ny, Nz, mpiRank);
   getRange(rhok, "rhok", Nx, Ny, Nz, mpiRank);
   getRange(rhok2, "rhok2", Nx, Ny, Nz, mpiRank);
   getRange(rho2, "rho2", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(rhoU1, "rhoU1", Nx, Ny, Nz, mpiRank);
   getRange(rhoUk, "rhoUk", Nx, Ny, Nz, mpiRank);
   getRange(rhoUk2, "rhoUk2", Nx, Ny, Nz, mpiRank);
   getRange(rhoU2, "rhoU2", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(rhoV1, "rhoV1", Nx, Ny, Nz, mpiRank);
   getRange(rhoVk, "rhoVk", Nx, Ny, Nz, mpiRank);
   getRange(rhoVk2, "rhoVk2", Nx, Ny, Nz, mpiRank);
   getRange(rhoV2, "rhoV2", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(rhoW1, "rhoW1", Nx, Ny, Nz, mpiRank);
   getRange(rhoWk, "rhoWk", Nx, Ny, Nz, mpiRank);
   getRange(rhoWk2, "rhoWk2", Nx, Ny, Nz, mpiRank);
   getRange(rhoW2, "rhoW2", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(rhoE1, "rhoE1", Nx, Ny, Nz, mpiRank);
   getRange(rhoEk, "rhoEk", Nx, Ny, Nz, mpiRank);
   getRange(rhoEk2, "rhoEk2", Nx, Ny, Nz, mpiRank);
   getRange(rhoE2, "rhoE2", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(p, "p", Nx, Ny, Nz, mpiRank);
   getRange(U, "U", Nx, Ny, Nz, mpiRank);
   getRange(V, "V", Nx, Ny, Nz, mpiRank);
   getRange(W, "W", Nx, Ny, Nz, mpiRank);
   getRange(T, "T", Nx, Ny, Nz, mpiRank);
   getRange(mu, "mu", Nx, Ny, Nz, mpiRank);
   getRange(sos, "sos", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;

}

