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

    c2d->allocY(rankfield);

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
    c2d->allocY(mom1_1);
    c2d->allocY(mom1_2);
    c2d->allocY(mom1_3);

    //27
    c2d->allocY(mom2_1);
    c2d->allocY(mom2_2);
    c2d->allocY(mom2_3);

    //30
    c2d->allocY(mom3_1);
    c2d->allocY(mom3_2);
    c2d->allocY(mom3_3);

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
    c2d->allocY(U);
    c2d->allocY(V);
    c2d->allocY(W);
    c2d->allocY(T);
    c2d->allocY(Ucurv);
    c2d->allocY(Vcurv);
    c2d->allocY(Wcurv);

    //72
    c2d->allocY(p);
    c2d->allocY(mu);
    c2d->allocY(sos);

    //84
    c2d->allocX(tempX1); tempXVec.push_back(tempX1);
    c2d->allocX(tempX2); tempXVec.push_back(tempX2);
    c2d->allocX(tempX3); tempXVec.push_back(tempX3);
    c2d->allocX(tempX4); tempXVec.push_back(tempX4);
    c2d->allocX(tempX5); tempXVec.push_back(tempX5);
    c2d->allocX(tempX6); tempXVec.push_back(tempX6);
    c2d->allocX(tempX7); tempXVec.push_back(tempX7);
    c2d->allocX(tempX8); tempXVec.push_back(tempX8);
    c2d->allocX(tempX9); tempXVec.push_back(tempX9);
    c2d->allocX(tempX10); tempXVec.push_back(tempX10);

    //94
    c2d->allocY(tempY1); tempYVec.push_back(tempY1);
    c2d->allocY(tempY2); tempYVec.push_back(tempY2);
    c2d->allocY(tempY3); tempYVec.push_back(tempY3);
    c2d->allocY(tempY4); tempYVec.push_back(tempY4);
    c2d->allocY(tempY5); tempYVec.push_back(tempY5);
    c2d->allocY(tempY6); tempYVec.push_back(tempY6);
    c2d->allocY(tempY7); tempYVec.push_back(tempY7);
    c2d->allocY(tempY8); tempYVec.push_back(tempY8);
    c2d->allocY(tempY9); tempYVec.push_back(tempY9);
    c2d->allocY(tempY10); tempYVec.push_back(tempY10);
    c2d->allocY(tempY11); tempYVec.push_back(tempY11);
    c2d->allocY(tempY12); tempYVec.push_back(tempY12);
    c2d->allocY(tempY13); tempYVec.push_back(tempY13);
    c2d->allocY(tempY14); tempYVec.push_back(tempY14);
    c2d->allocY(tempY15); tempYVec.push_back(tempY15);

    //104
    c2d->allocZ(tempZ1); tempZVec.push_back(tempZ1);
    c2d->allocZ(tempZ2); tempZVec.push_back(tempZ2);
    c2d->allocZ(tempZ3); tempZVec.push_back(tempZ3);
    c2d->allocZ(tempZ4); tempZVec.push_back(tempZ4);
    c2d->allocZ(tempZ5); tempZVec.push_back(tempZ5);
    c2d->allocZ(tempZ6); tempZVec.push_back(tempZ6);
    c2d->allocZ(tempZ7); tempZVec.push_back(tempZ7);
    c2d->allocZ(tempZ8); tempZVec.push_back(tempZ8);
    c2d->allocZ(tempZ9); tempZVec.push_back(tempZ9);
    c2d->allocZ(tempZ10); tempZVec.push_back(tempZ10);


    FOR_XYZ_YPEN{
	rankfield[ip] = mpiRank;
    }
    
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
	UChar_dx[ip] = J11[ip]*(fabs(U[ip]) + sos[ip])/dom->dx + J22[ip]*(fabs(V[ip])+sos[ip])/dom->dy + J33[ip]*(fabs(W[ip]) + sos[ip])/dom->dz;
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

    if(ts->timeSteppingType==Options::CONST_CFL){
	ts->dt = ts->CFL/max_UChar_dx;
    }else if(ts->timeSteppingType==Options::CONST_DT){
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

void CurvilinearCSolver::computeGradient(vector<double*> vecIn, vector<double*>vecOut){

    if(vecIn.size() > 5){
	cout << "ERROR:computeGradient: NEED MORE TEMPORARY MEMORY, COMPUTING GRADIENT ON LIST LARGER THAN 5" << endl;
    }else{

	int vecInSize = vecIn.size();

	if(vecOut.size() != vecInSize*3){
	    cout << "ERROR:computeGradient: vecOut size should be 3x the vecIn size!" << endl;
	}

	if(compStyle == VANILLA){

	    //Xi2 Derivatives First
	    for(int ip = 0; ip < vecInSize; ip++){
	        derivXi2->calc1stDerivField(vecIn[ip], vecOut[3*ip+1]);
       	    }

	    //Xi1 Derivatives Next
	    for(int ip = 0; ip < vecInSize; ip++){
	        c2d->transposeY2X_MajorIndex(vecIn[ip], tempXVec[ip]);
	    }

	    for(int ip = 0; ip < vecInSize; ip++){
	        derivXi1->calc1stDerivField(tempXVec[ip], tempXVec[ip+vecInSize]);
	    }

	    for(int ip = 0; ip < vecInSize; ip++){
	        c2d->transposeX2Y_MajorIndex(tempXVec[ip+vecInSize], vecOut[3*ip]);
	    }


	    //Xi3 Derivatives Finally
	    for(int ip = 0; ip < vecInSize; ip++){
	        c2d->transposeY2Z_MajorIndex(vecIn[ip], tempZVec[ip]);
	    }

	    for(int ip = 0; ip < vecInSize; ip++){
	        derivXi3->calc1stDerivField(tempZVec[ip], tempZVec[ip+vecInSize]);
	    }

	    for(int ip = 0; ip < vecInSize; ip++){
	        c2d->transposeZ2Y_MajorIndex(tempZVec[ip+vecInSize], vecOut[3*ip+2]);
	    }
	}else if(compStyle == OCC){

	    if(sbufVec.size() < 2*vecInSize){
		//Need to allocate memory in send buffer arrays...
		for(int ip = sbufVec.size(); ip < 2*vecInSize; ip++){
		    double *tmp_ptr = new double[c2d->decompBufSize];
		    sbufVec.push_back(tmp_ptr);
		}
	    }

	    if(rbufVec.size() < 2*vecInSize){
		//Need to allocate memory in recv buffer arrays...
		for(int ip = rbufVec.size(); ip < 2*vecInSize; ip++){
		    double *tmp_ptr = new double[c2d->decompBufSize];
		    rbufVec.push_back(tmp_ptr);
		}
	    }
	
	    //Going to need some handles for each of the communications
	    vector<MPI_Request> y2xHandles(vecInSize);
	    vector<MPI_Request> y2zHandles(vecInSize);

	    //Going to do all of our sends to x and z first
	    for(int ip = 0; ip < vecInSize; ip++){
		c2d->transposeY2X_MajorIndex_Start(y2xHandles[ip], vecIn[ip], tempXVec[ip], sbufVec[ip], rbufVec[ip]);
	    }

	    for(int ip = 0; ip < vecInSize; ip++){
		c2d->transposeY2Z_MajorIndex_Start(y2zHandles[ip], vecIn[ip], tempZVec[ip], sbufVec[vecInSize+ip], rbufVec[vecInSize+ip]);
	    }

	    //Now do our Xi2 Derivatives...
	    for(int ip = 0; ip < vecInSize; ip++){
	        derivXi2->calc1stDerivField(vecIn[ip], vecOut[3*ip+1]);
       	    }

     	    //Now lets wait for our first sends from Y2X, do our calculations and then send them back out
	    for(int ip = 0; ip < vecInSize; ip++){
		c2d->transposeY2X_MajorIndex_Wait(y2xHandles[ip], vecIn[ip], tempXVec[ip], sbufVec[ip], rbufVec[ip]);
		derivXi1->calc1stDerivField(tempXVec[ip], tempXVec[ip+vecInSize]);
		c2d->transposeX2Y_MajorIndex_Start(y2xHandles[ip], tempXVec[ip+vecInSize], vecOut[3*ip], sbufVec[ip], rbufVec[ip]);
	    }

	    //Now wait for our second sends from Y2Z, do our calculations and then send them back out 
	    for(int ip = 0; ip < vecInSize; ip++){
		c2d->transposeY2Z_MajorIndex_Wait(y2zHandles[ip], vecIn[ip], tempZVec[ip], sbufVec[vecInSize+ip], rbufVec[vecInSize+ip]);
		derivXi3->calc1stDerivField(tempZVec[ip], tempZVec[ip+vecInSize]);
		c2d->transposeZ2Y_MajorIndex_Start(y2zHandles[ip], tempZVec[ip+vecInSize], vecOut[3*ip+2], sbufVec[vecInSize+ip], rbufVec[vecInSize+ip]);
	    }

	    //Collect the first transpose 
	    for(int ip = 0; ip < vecInSize; ip++){
		c2d->transposeX2Y_MajorIndex_Wait(y2xHandles[ip], tempXVec[ip+vecInSize], vecOut[3*ip], sbufVec[ip], rbufVec[ip]);
	    }

	    //Collect the second transpose
	    for(int ip = 0; ip < vecInSize; ip++){
		c2d->transposeZ2Y_MajorIndex_Wait(y2zHandles[ip], tempZVec[ip+vecInSize], vecOut[3*ip+2], sbufVec[vecInSize+ip], rbufVec[vecInSize+ip]);
	    }

	}

    }

};

void CurvilinearCSolver::computeGradDotComponents(vector<double*> vecIn, vector<double*>vecOut){

    if(vecIn.size()/3 > 5){
	cout << "ERROR:computeGradDotComponents: NEED MORE TEMPORARY MEMORY, COMPUTING GRADIENT ON LIST LARGER THAN 5" << endl;
    }else{

	int vecInSize = vecIn.size();
	int N = vecIn.size()/3;

	if(vecOut.size() != vecInSize){
	    cout << "ERROR:computeGradDotComponents: vecOut size should be the same as the vecIn size!" << endl;
	}

	if(N*3 != vecInSize){
	    cout << "ERROR::computeGradDotComponents: missing an awkward number of inputs/outputs" << endl;
	}

        if(compStyle == VANILLA){

	    //Xi2 Derivatives First
	    for(int ip = 0; ip < N; ip++){
	        derivXi2->calc1stDerivField(vecIn[ip+N], vecOut[ip+N]);
       	    }

	    //Xi1 Derivatives Next
	    for(int ip = 0; ip < N; ip++){
	        c2d->transposeY2X_MajorIndex(vecIn[ip], tempXVec[ip]);
	    }

	    for(int ip = 0; ip < N; ip++){
	        derivXi1->calc1stDerivField(tempXVec[ip], tempXVec[ip+N]);
	    }

	    for(int ip = 0; ip < N; ip++){
	        c2d->transposeX2Y_MajorIndex(tempXVec[ip+N], vecOut[ip]);
	    }


	    //Xi3 Derivatives Finally
	    for(int ip = 0; ip < N; ip++){
	        c2d->transposeY2Z_MajorIndex(vecIn[2*N+ip], tempZVec[ip]);
	    }

	    for(int ip = 0; ip < N; ip++){
	        derivXi3->calc1stDerivField(tempZVec[ip], tempZVec[ip+N]);
	    }

	    for(int ip = 0; ip < N; ip++){
	        c2d->transposeZ2Y_MajorIndex(tempZVec[ip+N], vecOut[2*N+ip]);
	    }
	}else if(compStyle == OCC){

	
	    if(sbufVec.size() < 2*N){
		//Need to allocate memory in send buffer arrays...
		for(int ip = sbufVec.size(); ip < 2*N; ip++){
		    double *tmp_ptr = new double[c2d->decompBufSize];
		    sbufVec.push_back(tmp_ptr);
		}
	    }

	    if(rbufVec.size() < 2*N){
		//Need to allocate memory in recv buffer arrays...
		for(int ip = rbufVec.size(); ip < 2*N; ip++){
		    double *tmp_ptr = new double[c2d->decompBufSize];
		    rbufVec.push_back(tmp_ptr);
		}
	    }
	
	    //Going to need some handles for each of the communications
	    vector<MPI_Request> y2xHandles(vecInSize);
	    vector<MPI_Request> y2zHandles(vecInSize);

	    //Do our sends to x and z first
	    for(int ip = 0; ip < N; ip++){
		c2d->transposeY2X_MajorIndex_Start(y2xHandles[ip], vecIn[ip], tempXVec[ip], sbufVec[ip], rbufVec[ip]);
	    }

	    for(int ip = 0; ip < N; ip++){
		c2d->transposeY2Z_MajorIndex_Start(y2zHandles[ip], vecIn[2*N+ip], tempZVec[ip], sbufVec[ip+N], rbufVec[ip+N]);
	    }	

	   //Do our Xi2 Derivative Calculations
	   for(int ip = 0; ip < N; ip++){
	       derivXi2->calc1stDerivField(vecIn[ip+N], vecOut[ip+N]);
	   }

	   //Now we'll collect our transpose to x and do the calcs and send it back out
	    for(int ip = 0; ip < N; ip++){
		c2d->transposeY2X_MajorIndex_Wait(y2xHandles[ip], vecIn[ip], tempXVec[ip], sbufVec[ip], rbufVec[ip]);
		derivXi1->calc1stDerivField(tempXVec[ip], tempXVec[ip+N]);
		c2d->transposeX2Y_MajorIndex_Start(y2xHandles[ip], tempXVec[ip+N], vecOut[ip], sbufVec[ip], rbufVec[ip]);
	    }
	
	    //Do the same with out transposes to the Z
	    for(int ip = 0; ip < N; ip++){
		c2d->transposeY2Z_MajorIndex_Wait(y2zHandles[ip], vecIn[2*N+ip], tempZVec[ip], sbufVec[ip+N], rbufVec[ip+N]);
	        derivXi3->calc1stDerivField(tempZVec[ip], tempZVec[ip+N]);
		c2d->transposeZ2Y_MajorIndex_Start(y2zHandles[ip], tempZVec[ip+N], vecOut[2*N+ip], sbufVec[ip+N], rbufVec[ip+N]);

	    }

	    //Collect the final x->y transposes...
	    for(int ip = 0; ip < N; ip++){
		c2d->transposeX2Y_MajorIndex_Wait(y2xHandles[ip], tempXVec[ip+N], vecOut[ip], sbufVec[ip], rbufVec[ip]);
	    }

	    //Collect the final z->y transposes...
	    for(int ip = 0; ip < N; ip++){
		c2d->transposeZ2Y_MajorIndex_Wait(y2zHandles[ip], tempZVec[ip+N], vecOut[2*N+ip], sbufVec[ip+N], rbufVec[ip+N]);
	    }

	}

    }

};

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

    /////////////////////////////
    // Calculate the Gradients //
    /////////////////////////////

    double *vi[] = {U, V, W, T};
    vector<double*> vecIn(vi, vi+sizeof(vi)/sizeof(vi[0]));
    
    double *vo[] = {dU1, dU2, dU3, dV1, dV2, dV3, dW1, dW2, dW3, dT1, dT2, dT3};
    vector<double*> vecOut(vo, vo+sizeof(vo)/sizeof(vo[0]));

    computeGradient(vecIn, vecOut);

    if(useTiming){
	ftt2 = MPI_Wtime();
	IF_RANK0 cout << " > xidervtrans Timing: " << setw(6)  << (int)((ftt2-ftt1)*1000) << "ms" << endl;
	ftt1 = MPI_Wtime();
    }

    /////////////////////////////////////////
    // CALC TAU COMPONENTS & EQN COMPONENTS//
    /////////////////////////////////////////

    double *preCont_1, *preCont_2, *preCont_3;
    double *preMom1_1, *preMom1_2, *preMom1_3; 
    double *preMom2_1, *preMom2_2, *preMom2_3; 
    double *preMom3_1, *preMom3_2, *preMom3_3; 
    double *preEngy_1, *preEngy_2, *preEngy_3; 

    preMom1_1 = tempY1;  preMom1_2 = tempY2;  preMom1_3 = tempY3;
    preMom2_1 = tempY4;  preMom2_2 = tempY5;  preMom2_3 = tempY6;
    preMom3_1 = tempY7;  preMom3_2 = tempY8;  preMom3_3 = tempY9;
    preEngy_1 = tempY10; preEngy_2 = tempY11; preEngy_3 = tempY12;
    preCont_1 = tempY13; preCont_2 = tempY14; preCont_3 = tempY15;
    

    double b_1, b_2, b_3;
    double q_1, q_2, q_3;

    //Now recalculate properties in the new space
    FOR_XYZ_YPEN{

	//Viscous stress tensor
	Tau11[ip] =  (4.0/3.0)*(J11[ip]*dU1[ip] + J21[ip]*dU2[ip] + J31[ip]*dU3[ip]) + \
		    -(2.0/3.0)*(J12[ip]*dV1[ip] + J22[ip]*dV2[ip] + J32[ip]*dV3[ip]) + \
		    -(2.0/3.0)*(J13[ip]*dW1[ip] + J23[ip]*dW2[ip] + J33[ip]*dW3[ip]);

	Tau22[ip] = -(2.0/3.0)*(J11[ip]*dU1[ip] + J21[ip]*dU2[ip] + J31[ip]*dU3[ip]) + \
		     (4.0/3.0)*(J12[ip]*dV1[ip] + J22[ip]*dV2[ip] + J32[ip]*dV3[ip]) + \
		    -(2.0/3.0)*(J13[ip]*dW1[ip] + J23[ip]*dW2[ip] + J33[ip]*dW3[ip]);

	Tau33[ip] = -(2.0/3.0)*(J11[ip]*dU1[ip] + J21[ip]*dU2[ip] + J31[ip]*dU3[ip]) + \
		    -(2.0/3.0)*(J12[ip]*dV1[ip] + J22[ip]*dV2[ip] + J32[ip]*dV3[ip]) + \
		     (4.0/3.0)*(J13[ip]*dW1[ip] + J23[ip]*dW2[ip] + J33[ip]*dW3[ip]);

	Tau12[ip] =  J12[ip]*dU1[ip] + J22[ip]*dU2[ip] + J32[ip]*dU3[ip] + \
		     J11[ip]*dV1[ip] + J21[ip]*dV2[ip] + J31[ip]*dV3[ip];

	Tau13[ip] =  J13[ip]*dU1[ip] + J23[ip]*dU2[ip] + J33[ip]*dU3[ip] + \
		     J11[ip]*dW1[ip] + J21[ip]*dW2[ip] + J31[ip]*dW3[ip];

	Tau23[ip] =  J13[ip]*dV1[ip] + J23[ip]*dV2[ip] + J33[ip]*dV3[ip] + \
		     J12[ip]*dW1[ip] + J22[ip]*dW2[ip] + J32[ip]*dW3[ip];

	Tau11[ip] *= mu[ip];
	Tau22[ip] *= mu[ip];
	Tau33[ip] *= mu[ip];
	Tau12[ip] *= mu[ip];
	Tau13[ip] *= mu[ip];
	Tau23[ip] *= mu[ip];

	//Thermal conduction terms
	q_1 = (ig->cp/ig->Pr)*mu[ip]*(J11[ip]*dT1[ip] + J21[ip]*dT2[ip] + J31[ip]*dT3[ip]);
	q_2 = (ig->cp/ig->Pr)*mu[ip]*(J12[ip]*dT1[ip] + J22[ip]*dT2[ip] + J32[ip]*dT3[ip]);
	q_3 = (ig->cp/ig->Pr)*mu[ip]*(J13[ip]*dT1[ip] + J23[ip]*dT2[ip] + J33[ip]*dT3[ip]);

	//Viscous work + thermal conduction terms
	b_1 = U[ip]*Tau11[ip] + V[ip]*Tau12[ip] + W[ip]*Tau13[ip] + q_1;
	b_2 = U[ip]*Tau21[ip] + V[ip]*Tau22[ip] + W[ip]*Tau23[ip] + q_2;
	b_3 = U[ip]*Tau31[ip] + V[ip]*Tau32[ip] + W[ip]*Tau33[ip] + q_3;


	//Continuity Terms
	preCont_1[ip] = (1.0/J[ip])*(rhoP[ip]*Ucurv[ip]);
	preCont_2[ip] = (1.0/J[ip])*(rhoP[ip]*Vcurv[ip]);
	preCont_3[ip] = (1.0/J[ip])*(rhoP[ip]*Wcurv[ip]);

	//Combined Euler + viscous terms
	double F1, G1;
	double F2, G2;
	double F3, G3;

	F1 = rhoUP[ip]*Ucurv[ip] + J11[ip]*p[ip];
	F2 = rhoUP[ip]*Vcurv[ip] + J21[ip]*p[ip];
	F3 = rhoUP[ip]*Wcurv[ip] + J31[ip]*p[ip];

	G1 = J11[ip]*Tau11[ip] + J12[ip]*Tau21[ip] + J13[ip]*Tau31[ip];
	G2 = J21[ip]*Tau11[ip] + J22[ip]*Tau21[ip] + J23[ip]*Tau31[ip];
	G3 = J31[ip]*Tau11[ip] + J32[ip]*Tau21[ip] + J33[ip]*Tau31[ip];

	preMom1_1[ip] = (1.0/J[ip])*(F1 - G1);
	preMom1_2[ip] = (1.0/J[ip])*(F2 - G2);
	preMom1_3[ip] = (1.0/J[ip])*(F3 - G3);



	F1 = rhoVP[ip]*Ucurv[ip] + J12[ip]*p[ip];
	F2 = rhoVP[ip]*Vcurv[ip] + J22[ip]*p[ip];
	F3 = rhoVP[ip]*Wcurv[ip] + J32[ip]*p[ip];

	G1 = J11[ip]*Tau12[ip] + J12[ip]*Tau22[ip] + J13[ip]*Tau32[ip];
	G2 = J21[ip]*Tau12[ip] + J22[ip]*Tau22[ip] + J23[ip]*Tau32[ip];
	G3 = J31[ip]*Tau12[ip] + J32[ip]*Tau22[ip] + J33[ip]*Tau32[ip];

	preMom2_1[ip] = (1.0/J[ip])*(F1 - G1);
	preMom2_2[ip] = (1.0/J[ip])*(F2 - G2);
	preMom2_3[ip] = (1.0/J[ip])*(F3 - G3);



	F1 = rhoWP[ip]*Ucurv[ip] + J13[ip]*p[ip];
	F2 = rhoWP[ip]*Vcurv[ip] + J23[ip]*p[ip];
	F3 = rhoWP[ip]*Wcurv[ip] + J33[ip]*p[ip];

	G1 = J11[ip]*Tau13[ip] + J12[ip]*Tau23[ip] + J13[ip]*Tau33[ip];
	G2 = J21[ip]*Tau13[ip] + J22[ip]*Tau23[ip] + J23[ip]*Tau33[ip];
	G3 = J31[ip]*Tau13[ip] + J32[ip]*Tau23[ip] + J33[ip]*Tau33[ip];

	preMom3_1[ip] = (1.0/J[ip])*(F1 - G1);
	preMom3_2[ip] = (1.0/J[ip])*(F2 - G2);
	preMom3_3[ip] = (1.0/J[ip])*(F3 - G3);



	F1 = (rhoEP[ip] + p[ip])*Ucurv[ip];
	F2 = (rhoEP[ip] + p[ip])*Vcurv[ip];
	F3 = (rhoEP[ip] + p[ip])*Wcurv[ip];

	G1 = J11[ip]*b_1 + J12[ip]*b_2 + J13[ip]*b_3;
	G2 = J21[ip]*b_1 + J22[ip]*b_2 + J23[ip]*b_3;
	G3 = J31[ip]*b_1 + J32[ip]*b_2 + J33[ip]*b_3;


	preEngy_1[ip] = (1.0/J[ip])*(F1 - G1);
	preEngy_2[ip] = (1.0/J[ip])*(F2 - G2);
	preEngy_3[ip] = (1.0/J[ip])*(F3 - G3);

    }

    if(useTiming){
	ftt2 = MPI_Wtime();
	IF_RANK0 cout << " > Tau&Components Timing: " << setw(6)  << (int)((ftt2-ftt1)*1000) << "ms" << endl;
    }

    ///////////////////////////////////////
    // COMPUTE DERIVATIVES OF COMPONENTS //
    ///////////////////////////////////////

    double *vi2[] = {preCont_1, preMom1_1, preMom2_1, preMom3_1, preEngy_1, \
		     preCont_2, preMom1_2, preMom2_2, preMom3_2, preEngy_2, \
		     preCont_3, preMom1_3, preMom2_3, preMom3_3, preEngy_3};
    vector<double*> vecIn2(vi2, vi2+sizeof(vi2)/sizeof(vi2[0]));
    
    double *vo2[] = {cont_1, mom1_1, mom2_1, mom3_1, engy_1, \
		     cont_2, mom1_2, mom2_2, mom3_2, engy_2, \
		     cont_3, mom1_3, mom2_3, mom3_3, engy_3};
    vector<double*> vecOut2(vo2, vo2+sizeof(vo2)/sizeof(vo2[0]));

    computeGradDotComponents(vecIn2, vecOut2);


    if(useTiming){
	ftt2 = MPI_Wtime();
	IF_RANK0 cout << " >  derivtrans Timing: " << setw(6)  << (int)((ftt2-ftt1)*1000) << "ms" << endl;
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
		
	rhok2[ip]  = ts->dt*J[ip]*(-cont_1[ip] - cont_2[ip] - cont_3[ip] + (spgSource+contRHSSource(ip))/J[ip]);
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

	rhoUk2[ip] += -mom1_1[ip] - mom1_2[ip] -mom1_3[ip] + (spgSource+xmomRHSSource(ip))/J[ip];
	rhoUk2[ip] *= ts->dt*J[ip];
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

	rhoVk2[ip] += -mom2_1[ip] -mom2_2[ip] -mom2_3[ip] + (spgSource+ymomRHSSource(ip))/J[ip];
        rhoVk2[ip] *= ts->dt*J[ip];
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

	rhoWk2[ip] += -mom3_1[ip] -mom3_2[ip] -mom3_3[ip] + (spgSource+zmomRHSSource(ip))/J[ip];
        rhoWk2[ip] *= ts->dt*J[ip];
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

	rhoEk2[ip] = -engy_1[ip] - engy_2[ip] - engy_3[ip] + (spgSource+engyRHSSource(ip))/J[ip];
	rhoEk2[ip] *= ts->dt*J[ip];
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
            Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip];
            Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip];
            Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip];

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
            Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip];
            Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip];
            Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip];
	}
    }


    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > updNonCons Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
/*
    cout << "RHO" << endl; 
    FOR_Z_YPEN{
        FOR_Y_YPEN{
            FOR_X_YPEN{
                int ip = GETMAJIND_YPEN;
	 	cout << rhok[ip] << " ";
            }
	    cout << endl;
        }
	cout << endl;
    } 


    cout << "RHOU" << endl; 
    FOR_Z_YPEN{
        FOR_Y_YPEN{
            FOR_X_YPEN{
                int ip = GETMAJIND_YPEN;
	 	cout << rhoUk[ip] << " ";
            }
	    cout << endl;
        }
	cout << endl;
    }
 
    cout << "RHOV" << endl; 
    FOR_Z_YPEN{
        FOR_Y_YPEN{
            FOR_X_YPEN{
                int ip = GETMAJIND_YPEN;
	 	cout << rhoVk[ip] << " ";
            }
	    cout << endl;
        }
	cout << endl;
    } 
*/
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

	if(spongeFlag)
	    getRange(spg->sigma, "SIGMA", Nx, Ny, Nz, mpiRank);

        getRange(rho1, "RHO", Nx, Ny, Nz, mpiRank);
        getRange(U, "U", Nx, Ny, Nz, mpiRank);
        getRange(V, "V", Nx, Ny, Nz, mpiRank);
        getRange(W, "W", Nx, Ny, Nz, mpiRank);
        getRange(Ucurv, "Ucurv", Nx, Ny, Nz, mpiRank);
        getRange(Vcurv, "Vcurv", Nx, Ny, Nz, mpiRank);
        getRange(Wcurv, "Wcurv", Nx, Ny, Nz, mpiRank);
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

	//Converting over to X-Major Matrices
	FOR_Z_YPEN{
	    FOR_Y_YPEN{
		FOR_X_YPEN{
		    int ip = GETMAJIND_YPEN;
		    int jp = GETIND_YPEN;

		    tempY1[jp] = rho1[ip];
		    tempY2[jp] = rhoU1[ip];
		    tempY3[jp] = rhoV1[ip];
		    tempY4[jp] = rhoW1[ip];
		    tempY5[jp] = rhoE1[ip];

		    tempY6[jp] = msh->x[ip];
		    tempY7[jp] = msh->y[ip];
		    tempY8[jp] = msh->z[ip];

		}
	    }
	}
	
	


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

	double cNx, cNy, cNz;
	cNx = (double)Nx;	
	cNy = (double)Ny;	
	cNz = (double)Nz;	

	c2d->writeScalar(fh, disp, 1, &cNx); 
	c2d->writeScalar(fh, disp, 1, &cNy); 
	c2d->writeScalar(fh, disp, 1, &cNz); 
	c2d->writeVar(fh, disp, baseDirection, tempY6);
	c2d->writeVar(fh, disp, baseDirection, tempY7);
	c2d->writeVar(fh, disp, baseDirection, tempY8);
	c2d->writeVar(fh, disp, baseDirection, tempY1);
	c2d->writeVar(fh, disp, baseDirection, tempY2);
	c2d->writeVar(fh, disp, baseDirection, tempY3);
	c2d->writeVar(fh, disp, baseDirection, tempY4);
	c2d->writeVar(fh, disp, baseDirection, tempY5);

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

void CurvilinearCSolver::addImageOutput(PngWriter *png){
    imageList.push_back(png);
}

void CurvilinearCSolver::writeImages(){

    if(useTiming) ft1 = MPI_Wtime();

    //Get the timestep string just in case...
    ostringstream timeStepStringT;
    timeStepStringT << timeStep;
    int zeroPad = 6;
    string timeStepString = string(zeroPad - (timeStepStringT.str()).length(), '0') + timeStepStringT.str();

    //iterating through the image list
    for(list<PngWriter*>::iterator iter=imageList.begin(); iter != imageList.end(); ++iter){

	(*iter)->timeStepString = timeStepString;

	if(timeStep == 0 && (*iter)->dumpInterval == 0){

	    //if we haven't gotten an interpolator for this pngWriter yet, lets generate it
	    if((*iter)->interpolatorFlag == false){
		//Call the function to generate it...
		generateImagePlane(*iter); 	
		(*iter)->interpolatorFlag = true;
	    }

	    //Do the write...
	    writePlaneImageForVariable(*iter);	    
	}


	if(timeStep%(*iter)->dumpInterval==0 && (*iter)->dumpInterval > 0){

	    //if we haven't gotten an interpolator for this pngWriter yet, lets generate it
	    if((*iter)->interpolatorFlag == false){
		//Call the function to generate it...
		generateImagePlane(*iter); 	
		(*iter)->interpolatorFlag = true;
	    }

	    //Do the write
	    writePlaneImageForVariable(*iter);	    
	}
    }

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > writeImg Timing: " << setw(6)  << endl;
    }	 



}

void CurvilinearCSolver::generateImagePlane(PngWriter *pw){

    int N1 = pw->nx, N2 = pw->ny; 
    int plane    = pw->planeInd;
    int fraction = pw->fraction;

    double (*pointList)[3] = NULL;

    if(plane == 0){

	pointList = new double[N1*N2][3];

	//Get the plane location coordinate
	double x_plane = fraction*(msh->x_max[0]-msh->x_min[0]) + msh->x_min[0];

	//Get the coordinates for this case...
	//each pixel should represent a equal area box
	
	//See which case is the limiting box size
	double d1 = (msh->x_max[1]-msh->x_min[1])/((double)N1+1.0);
	double d2 = (msh->x_max[2]-msh->x_min[2])/((double)N2+1.0);

	//limited by the larger box...
	double dx = fmax(d1, d2);

	//Doing the pixed kind of as a "cell centered" location
	double base1 = msh->x_min[1] + dx/2;
	double base2 = msh->x_min[2] + dx/2;

	//Finally calculating the positions of the pixels
	for(int ip = 0; ip < N1; ip++){
	    for(int jp = 0; jp < N2; jp++){
		int ii = ip*N2 + jp;
		pointList[ii][0] = x_plane;
		pointList[ii][1] = base1 + dx*(double)ip;
		pointList[ii][2] = base2 + dx*(double)jp;
	    }
	}

    }else if(plane == 1){

	pointList = new double[N1*N2][3];

	//Get the plane location coordinate
	double y_plane = fraction*(msh->x_max[1] - msh->x_min[1]) + msh->x_min[1];

	//See which case is the limited box size
	double d1 = (msh->x_max[2] - msh->x_min[2])/((double)N1 + 1.0);
	double d2 = (msh->x_max[0] - msh->x_min[0])/((double)N2 + 1.0);

	double dx = fmax(d1, d2);

	//Get the base locations
	double base1 = msh->x_min[2] + dx/2;
	double base2 = msh->x_min[0] + dx/2;

	//Now calculate the positions of the pixels
	for(int ip = 0; ip < N1; ip++){
	    for(int jp = 0; jp < N2; jp++){
		int ii = ip*N2 + jp;
		pointList[ii][0] = base2 + dx*(double)jp;
		pointList[ii][1] = y_plane;
		pointList[ii][2] = base1 + dx*(double)ip;

	    }
	}

    }else if(plane == 2){
 
	pointList = new double[N1*N2][3];

	//Get the plane location coordinate
	double z_plane = fraction*(msh->x_max[2]-msh->x_min[2]) + msh->x_min[2];

	//See which case is the limited box size
	double d1 = (msh->x_max[0] - msh->x_min[0])/((double)N1 + 1.0);
	double d2 = (msh->x_max[1] - msh->x_min[1])/((double)N2 + 1.0);

	double dx = fmax(d1, d2);

	//Get the base locations
	double base1 = msh->x_min[0] + dx/2;
	double base2 = msh->x_min[1] + dx/2;

	//Now calculate the positions of the pixels
	for(int ip = 0; ip < N1; ip++){
	    for(int jp = 0; jp < N2; jp++){
	        int ii = ip*N2 + jp;
		pointList[ii][0] = base1 + dx*(double)ip;
		pointList[ii][1] = base2 + dx*(double)jp;	
		pointList[ii][2] = z_plane;
	    }
	}

   }else{
	cout << "Shouldn't be here in image writer?" << endl;
   }
    
    pw->ci = new CurvilinearInterpolator(this, pointList, N1*N2);

    delete[] pointList;

}

void CurvilinearCSolver::writePlaneImageForVariable(PngWriter *pw){

    int N1 = pw->nx;
    int N2 = pw->ny;
    double *var = pw->fieldPtr;

    double *ff_ci = new double[N1*N2];
    pw->ci->interpolateData(var, ff_ci);
  
    double *ff_local = new double[N1*N2];
    double *ff       = new double[N1*N2];
    for(int ip = 0; ip < N1*N2; ip++){
	ff_local[ip] = -1000000.0;
    }
    
    for(int ip = 0; ip < pw->ci->localPointFoundCount; ip++){
	ff_local[pw->ci->pointIndex[ip]] = ff_ci[ip];
    }

    delete[] ff_ci;

    MPI_Reduce(ff_local, ff, N1*N2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    delete[] ff_local;

    IF_RANK0{

	double dataMin =  100000000.0;
	double dataMax = -100000000.0;


	if(pw->valFlag == true){
	    dataMin = pw->valMin;
	    dataMax = pw->valMax;
	}else{
	    //get the max/min value in the plane...
	    for(int ip = 0; ip < N1*N2; ip++){
	        if(ff[ip] > -100000.0){
	            dataMin = fmin(dataMin, ff[ip]);
	            dataMax = fmax(dataMax, ff[ip]);
	        }
	    }
	}

	//Scale pixel value to local min and max
	int *r = new int[N1*N2];
	int *g = new int[N1*N2];
	int *b = new int[N1*N2];
	for(int ip = 0; ip < N1*N2; ip++){

	   double phitemp = (ff[ip] - dataMin)/(dataMax - dataMin); 

	   if((dataMax - dataMin) < 1E-6){
		phitemp = 0.0;
	   }

	   if(phitemp > 1.0){
	       phitemp = 1.0;
	   }

	   if(phitemp < 0.0){
	       phitemp = 0.0;
	   }

	   if(pw->cm == PngWriter::BWR){
	       pw->getPARAVIEW_BWR(phitemp, r[ip], g[ip], b[ip]);
	   }else if(pw->cm == PngWriter::RAINBOW){
	       pw->getRainbowColormap(phitemp, r[ip], g[ip], b[ip]);
	   }else if(pw->cm == PngWriter::GREYSCALE){ 
	       r[ip] = (int)(phitemp*255.0);
	       g[ip] = (int)(phitemp*255.0);
	       b[ip] = (int)(phitemp*255.0);
	   }	  

	}

	for(int jp = 0; jp < N2; jp++){
	    for(int ip = 0; ip < N1; ip++){
		int ii = jp*N1 + ip;
		//Grayscale image...
		if(ff[ii] > -100000.0){
		    pw->set(ip, jp, r[ii], g[ii], b[ii]);
		}else{
		    pw->set(ip, jp, 73, 175, 205);
		}
	    }
	}

	string imageName = pw->varName;
	imageName.append(".");
	imageName.append(pw->timeStepString);
	imageName.append(".png");
	pw->write(imageName.c_str());

	delete[] r;
	delete[] g;
	delete[] b;
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

   getRange(dU1, "dU1", Nx, Ny, Nz, mpiRank);
   getRange(dU2, "dU2", Nx, Ny, Nz, mpiRank);
   getRange(dU3, "dU3", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(dV1, "dV1", Nx, Ny, Nz, mpiRank);
   getRange(dV2, "dV2", Nx, Ny, Nz, mpiRank);
   getRange(dV3, "dV3", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(dW1, "dW1", Nx, Ny, Nz, mpiRank);
   getRange(dW2, "dW2", Nx, Ny, Nz, mpiRank);
   getRange(dW3, "dW3", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(dT1, "dT1", Nx, Ny, Nz, mpiRank);
   getRange(dT2, "dT2", Nx, Ny, Nz, mpiRank);
   getRange(dT3, "dT3", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(Tau11, "Tau11", Nx, Ny, Nz, mpiRank);
   getRange(Tau22, "Tau22", Nx, Ny, Nz, mpiRank);
   getRange(Tau33, "Tau33", Nx, Ny, Nz, mpiRank);
   getRange(Tau12, "Tau12", Nx, Ny, Nz, mpiRank);
   getRange(Tau23, "Tau23", Nx, Ny, Nz, mpiRank);
   getRange(Tau13, "Tau13", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(cont_1, "cont_1", Nx, Ny, Nz, mpiRank);
   getRange(cont_2, "cont_2", Nx, Ny, Nz, mpiRank);
   getRange(cont_3, "cont_3", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(mom1_1, "mom1_1", Nx, Ny, Nz, mpiRank);
   getRange(mom1_2, "mom1_2", Nx, Ny, Nz, mpiRank);
   getRange(mom1_3, "mom1_3", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(mom2_1, "mom2_1", Nx, Ny, Nz, mpiRank);
   getRange(mom2_2, "mom2_2", Nx, Ny, Nz, mpiRank);
   getRange(mom2_3, "mom2_3", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(mom3_1, "mom3_1", Nx, Ny, Nz, mpiRank);
   getRange(mom3_2, "mom3_2", Nx, Ny, Nz, mpiRank);
   getRange(mom3_3, "mom3_3", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;
   getRange(engy_1, "engy_1", Nx, Ny, Nz, mpiRank);
   getRange(engy_2, "engy_2", Nx, Ny, Nz, mpiRank);
   getRange(engy_3, "engy_3", Nx, Ny, Nz, mpiRank);
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
   getRange(Ucurv, "Ucurv", Nx, Ny, Nz, mpiRank);
   getRange(Vcurv, "Vcurv", Nx, Ny, Nz, mpiRank);
   getRange(Wcurv, "Wcurv", Nx, Ny, Nz, mpiRank);
   getRange(T, "T", Nx, Ny, Nz, mpiRank);
   getRange(mu, "mu", Nx, Ny, Nz, mpiRank);
   getRange(sos, "sos", Nx, Ny, Nz, mpiRank);
   IF_RANK0 cout << " " << endl;

}


