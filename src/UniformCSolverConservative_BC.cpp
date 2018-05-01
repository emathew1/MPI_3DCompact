#include "UniformCSolverConservative.hpp"

void UniformCSolverConservative::setInitialConditions(){

    if(useTiming) ft1 = MPI_Wtime();

    IF_RANK0{
        cout << endl;
        cout << " > Setting initial conditions..." << endl; 
    }

    //just do the simple stuff in a loop...
    FOR_XYZ_XPEN{
	rho1[ip] = rho0[ip];
	U[ip]	 = U0[ip];
	V[ip] 	 = V0[ip];
	W[ip] 	 = W0[ip];
	p[ip]	 = p0[ip];
	rhoU1[ip] = rho1[ip]*U[ip];	
	rhoV1[ip] = rho1[ip]*V[ip];	
	rhoW1[ip] = rho1[ip]*W[ip];
    }

    //Can now release initial condition data...
    c2d->deallocXYZ(rho0);
    c2d->deallocXYZ(U0);
    c2d->deallocXYZ(V0);
    c2d->deallocXYZ(W0);
    c2d->deallocXYZ(p0);

    //Call the ideal gas relations for the slightly more involved stuff..
    FOR_XYZ_XPEN{
	rhoE1[ip] = ig->solverhoE(rho1[ip], p[ip], U[ip], V[ip], W[ip]);
    	T[ip]     = ig->solveT(rho1[ip], p[ip]);
        mu[ip]    = ig->solveMu(T[ip]);
        Amu[ip]   = ig->solveAmu(T[ip]);
        sos[ip]   = ig->solveSOS(rho1[ip], p[ip]);
    }
    //This is where we'll do the boundary condition specific stuff...
    bool wallBCFlag = false;

    //--------------------------------
    //No-Slip Wall Boundary Conditions
    //--------------------------------

    if(bc->bcX0 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_X0_XPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
            T[ip] = calcNeumann(T[GETMAJIND_XPEN_Xp1],
                                T[GETMAJIND_XPEN_Xp2],
                                T[GETMAJIND_XPEN_Xp3],
                                T[GETMAJIND_XPEN_Xp4],
                                T[GETMAJIND_XPEN_Xp5],
                                T[GETMAJIND_XPEN_Xp6]);
        }END_FORX0
    }


    if(bc->bcX1 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_X1_XPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
            T[ip] = calcNeumann(T[GETMAJIND_XPEN_Xm1],
                                T[GETMAJIND_XPEN_Xm2],
                                T[GETMAJIND_XPEN_Xm3],
                                T[GETMAJIND_XPEN_Xm4],
                                T[GETMAJIND_XPEN_Xm5],
                                T[GETMAJIND_XPEN_Xm6]);
        }END_FORX1
    }

    if(bc->bcY0 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Y0_XPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
	    if(neumannLocalY){
	        T[ip] = calcNeumann(T[GETYPIND_XPEN_Yp1],
                                    T[GETYPIND_XPEN_Yp2],
                                    T[GETYPIND_XPEN_Yp3],
                                    T[GETYPIND_XPEN_Yp4],
                                    T[GETYPIND_XPEN_Yp5],
                                    T[GETYPIND_XPEN_Yp6]);
	    } 
        }END_FORY0

	if(!neumannLocalY){
	    double *T2;
	    c2d->allocY(T2);
	    c2d->transposeX2Y_MajorIndex(T, T2);

	    FOR_Y0_YPEN_MAJ{
                T2[ip] = calcNeumann(T2[GETMAJIND_YPEN_Yp1],
                                     T2[GETMAJIND_YPEN_Yp2],
                                     T2[GETMAJIND_YPEN_Yp3],
                                     T2[GETMAJIND_YPEN_Yp4],
                                     T2[GETMAJIND_YPEN_Yp5],
                                     T2[GETMAJIND_YPEN_Yp6]);

	    }END_FORY0 

	    c2d->transposeY2X_MajorIndex(T2, T);
	    c2d->deallocXYZ(T2);
	}
    }

    if(bc->bcY1 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Y1_XPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
	    if(neumannLocalY){
	        T[ip] = calcNeumann(T[GETYPIND_XPEN_Ym1],
                                    T[GETYPIND_XPEN_Ym2],
                                    T[GETYPIND_XPEN_Ym3],
                                    T[GETYPIND_XPEN_Ym4],
                                    T[GETYPIND_XPEN_Ym5],
                                    T[GETYPIND_XPEN_Ym6]);
	    } 
        }END_FORY1

	if(!neumannLocalY){
	    double *T2;
	    c2d->allocY(T2);
	    c2d->transposeX2Y_MajorIndex(T, T2);

	    FOR_Y1_YPEN_MAJ{
                T2[ip] = calcNeumann(T2[GETMAJIND_YPEN_Ym1],
                                     T2[GETMAJIND_YPEN_Ym2],
                                     T2[GETMAJIND_YPEN_Ym3],
                                     T2[GETMAJIND_YPEN_Ym4],
                                     T2[GETMAJIND_YPEN_Ym5],
                                     T2[GETMAJIND_YPEN_Ym6]);
 	    }END_FORY1 

	    c2d->transposeY2X_MajorIndex(T2, T);
	    c2d->deallocXYZ(T2);
	}
    }

    if(bc->bcZ0 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Z0_XPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
	    if(neumannLocalZ){
	        T[ip] = calcNeumann(T[GETZPIND_XPEN_Zp1],
                                    T[GETZPIND_XPEN_Zp2],
                                    T[GETZPIND_XPEN_Zp3],
                                    T[GETZPIND_XPEN_Zp4],
                                    T[GETZPIND_XPEN_Zp5],
                                    T[GETZPIND_XPEN_Zp6]);
	    } 
        }END_FORZ0

	if(!neumannLocalZ){
	    double *T2, *T3;
	    c2d->allocY(T2);
	    c2d->allocZ(T3);
	    c2d->transposeX2Y_MajorIndex(T, T2);
	    c2d->transposeY2Z_MajorIndex(T2, T3);

	    FOR_Z0_ZPEN_MAJ{
                T3[ip] = calcNeumann(T3[GETMAJIND_ZPEN_Zp1],
                                     T3[GETMAJIND_ZPEN_Zp2],
                                     T3[GETMAJIND_ZPEN_Zp3],
                                     T3[GETMAJIND_ZPEN_Zp4],
                                     T3[GETMAJIND_ZPEN_Zp5],
                                     T3[GETMAJIND_ZPEN_Zp6]);
	    }END_FORZ0

	    c2d->transposeZ2Y_MajorIndex(T3, T2);
	    c2d->transposeY2X_MajorIndex(T2, T);

	    c2d->deallocXYZ(T2);
	    c2d->deallocXYZ(T3);
	}
    }

    if(bc->bcZ1 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Z1_XPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
	    if(neumannLocalZ){
	        T[ip] = calcNeumann(T[GETZPIND_XPEN_Zm1],
                                    T[GETZPIND_XPEN_Zm2],
                                    T[GETZPIND_XPEN_Zm3],
                                    T[GETZPIND_XPEN_Zm4],
                                    T[GETZPIND_XPEN_Zm5],
                                    T[GETZPIND_XPEN_Zm6]);
	    } 

        }END_FORZ1

	if(!neumannLocalZ){
	    double *T2, *T3;
	    c2d->allocY(T2);
	    c2d->allocZ(T3);
	    c2d->transposeX2Y_MajorIndex(T, T2);
	    c2d->transposeY2Z_MajorIndex(T2, T3);

	    FOR_Z1_ZPEN_MAJ{
                T3[ip] = calcNeumann(T3[GETMAJIND_ZPEN_Zm1],
                                     T3[GETMAJIND_ZPEN_Zm2],
                                     T3[GETMAJIND_ZPEN_Zm3],
                                     T3[GETMAJIND_ZPEN_Zm4],
                                     T3[GETMAJIND_ZPEN_Zm5],
                                     T3[GETMAJIND_ZPEN_Zm6]);
 	    }END_FORZ1

	    c2d->transposeZ2Y_MajorIndex(T3, T2);
	    c2d->transposeY2X_MajorIndex(T2, T);

	    c2d->deallocXYZ(T2);
	    c2d->deallocXYZ(T3);
	}


    }

    //-------------------------------
    //Moving Wall Boundary Conditions
    //-------------------------------

    if(bc->bcX0 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_X0_XPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = X0WallV;
            W[ip]  = X0WallW;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = rho1[ip]*X0WallV;
            rhoW1[ip] = rho1[ip]*X0WallW;
            T[ip] = calcNeumann(T[GETMAJIND_XPEN_Xp1],
                                T[GETMAJIND_XPEN_Xp2],
                                T[GETMAJIND_XPEN_Xp3],
                                T[GETMAJIND_XPEN_Xp4],
                                T[GETMAJIND_XPEN_Xp5],
                                T[GETMAJIND_XPEN_Xp6]);
        }END_FORX0
    }

    if(bc->bcX1 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_X1_XPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = X1WallV;
            W[ip]  = X1WallW;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = rho1[ip]*X1WallV;
            rhoW1[ip] = rho1[ip]*X1WallW;
            T[ip] = calcNeumann(T[GETMAJIND_XPEN_Xm1],
                                T[GETMAJIND_XPEN_Xm2],
                                T[GETMAJIND_XPEN_Xm3],
                                T[GETMAJIND_XPEN_Xm4],
                                T[GETMAJIND_XPEN_Xm5],
                                T[GETMAJIND_XPEN_Xm6]);
        }END_FORX1
    }

    if(bc->bcY0 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Y0_XPEN_MAJ{
            U[ip]  = Y0WallU;
            V[ip]  = 0.0;
            W[ip]  = Y0WallW;
            rhoU1[ip] = rho1[ip]*Y0WallU;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = rho1[ip]*Y0WallW;
	    if(neumannLocalY){
	        T[ip] = calcNeumann(T[GETYPIND_XPEN_Yp1],
                                    T[GETYPIND_XPEN_Yp2],
                                    T[GETYPIND_XPEN_Yp3],
                                    T[GETYPIND_XPEN_Yp4],
                                    T[GETYPIND_XPEN_Yp5],
                                    T[GETYPIND_XPEN_Yp6]);
	    } 
        }END_FORY0

	if(!neumannLocalY){
 	    double *T2;
	    c2d->allocY(T2);
	    c2d->transposeX2Y_MajorIndex(T, T2);

	    FOR_Y0_YPEN_MAJ{
                T2[ip] = calcNeumann(T2[GETMAJIND_YPEN_Yp1],
                                     T2[GETMAJIND_YPEN_Yp2],
                                     T2[GETMAJIND_YPEN_Yp3],
                                     T2[GETMAJIND_YPEN_Yp4],
                                     T2[GETMAJIND_YPEN_Yp5],
                                     T2[GETMAJIND_YPEN_Yp6]);
	    }END_FORY0 

	    c2d->transposeY2X_MajorIndex(T2, T);
	    c2d->deallocXYZ(T2);
	}


   }

    if(bc->bcY1 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Y1_XPEN_MAJ{
            U[ip]  = Y1WallU;
            V[ip]  = 0.0;
            W[ip]  = Y1WallW;
            rhoU1[ip] = rho1[ip]*Y1WallU;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = rho1[ip]*Y1WallW;
	    if(neumannLocalY){
	        T[ip] = calcNeumann(T[GETYPIND_XPEN_Ym1],
                                    T[GETYPIND_XPEN_Ym2],
                                    T[GETYPIND_XPEN_Ym3],
                                    T[GETYPIND_XPEN_Ym4],
                                    T[GETYPIND_XPEN_Ym5],
                                    T[GETYPIND_XPEN_Ym6]);
	    } 
        }END_FORY1

	if(!neumannLocalY){
	    double *T2;
	    c2d->allocY(T2);
	    c2d->transposeX2Y_MajorIndex(T, T2);

	    FOR_Y1_YPEN_MAJ{
                T2[ip] = calcNeumann(T2[GETMAJIND_YPEN_Ym1],
                                     T2[GETMAJIND_YPEN_Ym2],
                                     T2[GETMAJIND_YPEN_Ym3],
                                     T2[GETMAJIND_YPEN_Ym4],
                                     T2[GETMAJIND_YPEN_Ym5],
                                     T2[GETMAJIND_YPEN_Ym6]);
 	    }END_FORY1 

	    c2d->transposeY2X_MajorIndex(T2, T);
	    c2d->deallocXYZ(T2);
	}

 
   }

    if(bc->bcZ0 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Z0_XPEN_MAJ{
            U[ip]  = Z0WallU;
            V[ip]  = Z0WallV;
            W[ip]  = 0.0;
            rhoU1[ip] = rho1[ip]*Z0WallU;
            rhoV1[ip] = rho1[ip]*Z0WallV;
            rhoW1[ip] = 0.0;
	    if(neumannLocalZ){
	        T[ip] = calcNeumann(T[GETZPIND_XPEN_Zp1],
                                    T[GETZPIND_XPEN_Zp2],
                                    T[GETZPIND_XPEN_Zp3],
                                    T[GETZPIND_XPEN_Zp4],
                                    T[GETZPIND_XPEN_Zp5],
                                    T[GETZPIND_XPEN_Zp6]);
	    } 
        }END_FORZ0


	if(!neumannLocalZ){
	    double *T2, *T3;
	    c2d->allocY(T2);
	    c2d->allocZ(T3);
	    c2d->transposeX2Y_MajorIndex(T, T2);
	    c2d->transposeY2Z_MajorIndex(T2, T3);

	    FOR_Z0_ZPEN_MAJ{
                T3[ip] = calcNeumann(T3[GETMAJIND_ZPEN_Zp1],
                                     T3[GETMAJIND_ZPEN_Zp2],
                                     T3[GETMAJIND_ZPEN_Zp3],
                                     T3[GETMAJIND_ZPEN_Zp4],
                                     T3[GETMAJIND_ZPEN_Zp5],
                                     T3[GETMAJIND_ZPEN_Zp6]);
	    }END_FORZ0

	    c2d->transposeZ2Y_MajorIndex(T3, T2);
	    c2d->transposeY2X_MajorIndex(T2, T);

	    c2d->deallocXYZ(T2);
	    c2d->deallocXYZ(T3);
	}


    }

    if(bc->bcZ1 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Z1_XPEN_MAJ{
            U[ip]  = Z1WallU;
            V[ip]  = Z1WallV;
            W[ip]  = 0.0;
            rhoU1[ip] = rho1[ip]*Z1WallU;
            rhoV1[ip] = rho1[ip]*Z1WallV;
            rhoW1[ip] = 0.0;
	    if(neumannLocalZ){
	        T[ip] = calcNeumann(T[GETZPIND_XPEN_Zm1],
                                    T[GETZPIND_XPEN_Zm2],
                                    T[GETZPIND_XPEN_Zm3],
                                    T[GETZPIND_XPEN_Zm4],
                                    T[GETZPIND_XPEN_Zm5],
                                    T[GETZPIND_XPEN_Zm6]);
	    } 

        }END_FORZ1


	if(!neumannLocalZ){
	    double *T2, *T3;
	    c2d->allocY(T2);
	    c2d->allocZ(T3);
	    c2d->transposeX2Y_MajorIndex(T, T2);
	    c2d->transposeY2Z_MajorIndex(T2, T3);

	    FOR_Z1_ZPEN_MAJ{
                T3[ip] = calcNeumann(T3[GETMAJIND_ZPEN_Zm1],
                                     T3[GETMAJIND_ZPEN_Zm2],
                                     T3[GETMAJIND_ZPEN_Zm3],
                                     T3[GETMAJIND_ZPEN_Zm4],
                                     T3[GETMAJIND_ZPEN_Zm5],
                                     T3[GETMAJIND_ZPEN_Zm6]);
 	    }END_FORZ1

	    c2d->transposeZ2Y_MajorIndex(T3, T2);
	    c2d->transposeY2X_MajorIndex(T2, T);

	    c2d->deallocXYZ(T2);
	    c2d->deallocXYZ(T3);
	}

    }




    if(wallBCFlag == true){
	//Need to update the pressure, sos, and rhoE fields at the boundaries with walls...
	FOR_XYZ_XPEN{
	    p[ip]     = ig->solvep_idealgas(rho1[ip], T[ip]);
	    sos[ip]   = ig->solveSOS(rho1[ip],p[ip]);
	    rhoE1[ip] = ig->solverhoE(rho1[ip], p[ip], U[ip], V[ip], W[ip]);
	}
    }

    if(spongeFlag == true){
	FOR_XYZ_XPEN{
	    spg->spongeRhoAvg[ip]  = rho1[ip];
	    spg->spongeRhoUAvg[ip] = rhoU1[ip];
	    spg->spongeRhoVAvg[ip] = rhoV1[ip];
	    spg->spongeRhoWAvg[ip] = rhoW1[ip];
	    spg->spongeRhoEAvg[ip] = rhoE1[ip];
 	}
    }

    IF_RANK0 std::cout << " > Finished initialization of flow field " << std::endl;

    getRange(rho1,  "rho0", pxSize[0], pxSize[1], pxSize[2], mpiRank);
    getRange(rhoU1, "rhoU0", pxSize[0], pxSize[1], pxSize[2], mpiRank);
    getRange(rhoV1, "rhoV0", pxSize[0], pxSize[1], pxSize[2], mpiRank);
    getRange(rhoW1, "rhoW0", pxSize[0], pxSize[1], pxSize[2], mpiRank);
    getRange(rhoE1, "rhoE0", pxSize[0], pxSize[1], pxSize[2], mpiRank);
    getRange(T, "T0", pxSize[0], pxSize[1], pxSize[2], mpiRank);
    getRange(p, "p0", pxSize[0], pxSize[1], pxSize[2], mpiRank);

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > setInitCond Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }


}

void UniformCSolverConservative::preStepBCHandling(){

    if(useTiming) ft1 = MPI_Wtime();

    double *rhoP, *rhoUP, *rhoVP, *rhoWP, *rhoEP;
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


    //--------------------------------
    //No-slip wall boundary conditions
    //--------------------------------

    if(bc->bcX0 == BC::ADIABATIC_WALL){
	FOR_X0_XPEN_MAJ{
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    T[ip] = calcNeumann(T[GETMAJIND_XPEN_Xp1],
				T[GETMAJIND_XPEN_Xp2],
				T[GETMAJIND_XPEN_Xp3],
				T[GETMAJIND_XPEN_Xp4],
				T[GETMAJIND_XPEN_Xp5],
				T[GETMAJIND_XPEN_Xp6]);
	   p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	   rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
	}END_FORX0

    }
    

    if(bc->bcX1 == BC::ADIABATIC_WALL){
	FOR_X1_XPEN_MAJ{
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    T[ip] = calcNeumann(T[GETMAJIND_XPEN_Xm1],
				T[GETMAJIND_XPEN_Xm2],
				T[GETMAJIND_XPEN_Xm3],
				T[GETMAJIND_XPEN_Xm4],
				T[GETMAJIND_XPEN_Xm5],
				T[GETMAJIND_XPEN_Xm6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
	}END_FORX1
    }   


    if(bc->bcY0 == BC::ADIABATIC_WALL){

	if(!neumannLocalY){
            double *T2;
	    c2d->allocY(T2);
	    c2d->transposeX2Y_MajorIndex(T, T2);

	    FOR_Y0_YPEN_MAJ{
                T2[ip] = calcNeumann(T2[GETMAJIND_YPEN_Yp1],
                                     T2[GETMAJIND_YPEN_Yp2],
                                     T2[GETMAJIND_YPEN_Yp3],
                                     T2[GETMAJIND_YPEN_Yp4],
                                     T2[GETMAJIND_YPEN_Yp5],
                                     T2[GETMAJIND_YPEN_Yp6]);
	    }END_FORY0

	    c2d->transposeY2X_MajorIndex(T2, T);
	    c2d->deallocXYZ(T2);
	}

	FOR_Y0_XPEN_MAJ{
	    if(neumannLocalY){
	        T[ip] = calcNeumann(T[GETYPIND_XPEN_Yp1],
                                    T[GETYPIND_XPEN_Yp2],
                                    T[GETYPIND_XPEN_Yp3],
                                    T[GETYPIND_XPEN_Yp4],
                                    T[GETYPIND_XPEN_Yp5],
                                    T[GETYPIND_XPEN_Yp6]);
	    } 
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]);
	}END_FORY0
    }
    

    if(bc->bcY1 == BC::ADIABATIC_WALL){

	if(!neumannLocalY){
	    double *T2;
	    c2d->allocY(T2);
	    c2d->transposeX2Y_MajorIndex(T, T2);

	    FOR_Y1_YPEN_MAJ{
                T2[ip] = calcNeumann(T2[GETMAJIND_YPEN_Ym1],
                                     T2[GETMAJIND_YPEN_Ym2],
                                     T2[GETMAJIND_YPEN_Ym3],
                                     T2[GETMAJIND_YPEN_Ym4],
                                     T2[GETMAJIND_YPEN_Ym5],
                                     T2[GETMAJIND_YPEN_Ym6]);
 	    }END_FORY1 

	    c2d->transposeY2X_MajorIndex(T2, T);
	    c2d->deallocXYZ(T2);
	}

	FOR_Y1_XPEN_MAJ{
	    if(neumannLocalY){
	        T[ip] = calcNeumann(T[GETYPIND_XPEN_Ym1],
                                    T[GETYPIND_XPEN_Ym2],
                                    T[GETYPIND_XPEN_Ym3],
                                    T[GETYPIND_XPEN_Ym4],
                                    T[GETYPIND_XPEN_Ym5],
                                    T[GETYPIND_XPEN_Ym6]);
	    } 
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
	}END_FORY1
    }



    if(bc->bcZ0 == BC::ADIABATIC_WALL){

	if(!neumannLocalZ){
	    double *T2, *T3;
	    c2d->allocY(T2);
	    c2d->allocZ(T3);
	    c2d->transposeX2Y_MajorIndex(T, T2);
	    c2d->transposeY2Z_MajorIndex(T2, T3);

	    FOR_Z0_ZPEN_MAJ{
                T3[ip] = calcNeumann(T3[GETMAJIND_ZPEN_Zp1],
                                     T3[GETMAJIND_ZPEN_Zp2],
                                     T3[GETMAJIND_ZPEN_Zp3],
                                     T3[GETMAJIND_ZPEN_Zp4],
                                     T3[GETMAJIND_ZPEN_Zp5],
                                     T3[GETMAJIND_ZPEN_Zp6]);
	    }END_FORZ0

	    c2d->transposeZ2Y_MajorIndex(T3, T2);
	    c2d->transposeY2X_MajorIndex(T2, T);

	    c2d->deallocXYZ(T2);
	    c2d->deallocXYZ(T3);
	}
	

	FOR_Z0_XPEN{
	    if(neumannLocalZ){
	        T[ip] = calcNeumann(T[GETZPIND_XPEN_Zp1],
                                    T[GETZPIND_XPEN_Zp2],
                                    T[GETZPIND_XPEN_Zp3],
                                    T[GETZPIND_XPEN_Zp4],
                                    T[GETZPIND_XPEN_Zp5],
                                    T[GETZPIND_XPEN_Zp6]);
	    } 
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
	}END_FORZ0
    }



    if(bc->bcZ1 == BC::ADIABATIC_WALL){

	if(!neumannLocalZ){
	    double *T2, *T3;
	    c2d->allocY(T2);
	    c2d->allocZ(T3);
	    c2d->transposeX2Y_MajorIndex(T, T2);
	    c2d->transposeY2Z_MajorIndex(T2, T3);

	    FOR_Z1_ZPEN_MAJ{
                T3[ip] = calcNeumann(T3[GETMAJIND_ZPEN_Zm1],
                                     T3[GETMAJIND_ZPEN_Zm2],
                                     T3[GETMAJIND_ZPEN_Zm3],
                                     T3[GETMAJIND_ZPEN_Zm4],
                                     T3[GETMAJIND_ZPEN_Zm5],
                                     T3[GETMAJIND_ZPEN_Zm6]);
 	    }END_FORZ1

	    c2d->transposeZ2Y_MajorIndex(T3, T2);
	    c2d->transposeY2X_MajorIndex(T2, T);

	    c2d->deallocXYZ(T2);
	    c2d->deallocXYZ(T3);
	}

	FOR_Z1_XPEN{
	    if(neumannLocalZ){
	        T[ip] = calcNeumann(T[GETZPIND_XPEN_Zm1],
                                    T[GETZPIND_XPEN_Zm2],
                                    T[GETZPIND_XPEN_Zm3],
                                    T[GETZPIND_XPEN_Zm4],
                                    T[GETZPIND_XPEN_Zm5],
                                    T[GETZPIND_XPEN_Zm6]);
	    } 
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
	}END_FORZ1
    }


    //-------------------------------
    //Moving wall boundary conditions
    //-------------------------------

    if(bc->bcX0 == BC::MOVING_ADIABATIC_WALL){
        FOR_X0_XPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = X0WallV;
            W[ip]  = X0WallW;
            rhoUP[ip] = 0.0;
            rhoVP[ip] = rhoP[ip]*X0WallV;
            rhoWP[ip] = rhoP[ip]*X0WallW;
            T[ip] = calcNeumann(T[GETMAJIND_XPEN_Xp1],
                                T[GETMAJIND_XPEN_Xp2],
                                T[GETMAJIND_XPEN_Xp3],
                                T[GETMAJIND_XPEN_Xp4],
                                T[GETMAJIND_XPEN_Xp5],
                                T[GETMAJIND_XPEN_Xp6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
        }END_FORX0
    }

    if(bc->bcX1 == BC::MOVING_ADIABATIC_WALL){
        FOR_X1_XPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = X1WallV;
            W[ip]  = X1WallW;
            rhoUP[ip] = 0.0;
            rhoVP[ip] = rhoP[ip]*X1WallV;
            rhoWP[ip] = rhoP[ip]*X1WallW;
            T[ip] = calcNeumann(T[GETMAJIND_XPEN_Xm1],
                                T[GETMAJIND_XPEN_Xm2],
                                T[GETMAJIND_XPEN_Xm3],
                                T[GETMAJIND_XPEN_Xm4],
                                T[GETMAJIND_XPEN_Xm5],
                                T[GETMAJIND_XPEN_Xm6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
        }END_FORX1
    }

    if(bc->bcY0 == BC::MOVING_ADIABATIC_WALL){

	if(!neumannLocalY){
 	    double *T2;
	    c2d->allocY(T2);
	    c2d->transposeX2Y_MajorIndex(T, T2);

	    FOR_Y0_YPEN_MAJ{
                T2[ip] = calcNeumann(T2[GETMAJIND_YPEN_Yp1],
                                     T2[GETMAJIND_YPEN_Yp2],
                                     T2[GETMAJIND_YPEN_Yp3],
                                     T2[GETMAJIND_YPEN_Yp4],
                                     T2[GETMAJIND_YPEN_Yp5],
                                     T2[GETMAJIND_YPEN_Yp6]);
	    }END_FORY0 

	    c2d->transposeY2X_MajorIndex(T2, T);
	    c2d->deallocXYZ(T2);
	}
	
	
        FOR_Y0_XPEN{
	    if(neumannLocalY){
	        T[ip] = calcNeumann(T[GETYPIND_XPEN_Yp1],
                                    T[GETYPIND_XPEN_Yp2],
                                    T[GETYPIND_XPEN_Yp3],
                                    T[GETYPIND_XPEN_Yp4],
                                    T[GETYPIND_XPEN_Yp5],
                                    T[GETYPIND_XPEN_Yp6]);
	    } 
            U[ip]  = Y0WallU;
            V[ip]  = 0.0;
            W[ip]  = Y0WallW;
            rhoUP[ip] = rhoP[ip]*Y0WallU;
            rhoVP[ip] = 0.0;
            rhoWP[ip] = rhoP[ip]*Y0WallW;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
        }END_FORY0
    }

    if(bc->bcY1 == BC::MOVING_ADIABATIC_WALL){

	if(!neumannLocalY){
	    double *T2;
	    c2d->allocY(T2);
	    c2d->transposeX2Y_MajorIndex(T, T2);

	    FOR_Y1_YPEN_MAJ{
                T2[ip] = calcNeumann(T2[GETMAJIND_YPEN_Ym1],
                                     T2[GETMAJIND_YPEN_Ym2],
                                     T2[GETMAJIND_YPEN_Ym3],
                                     T2[GETMAJIND_YPEN_Ym4],
                                     T2[GETMAJIND_YPEN_Ym5],
                                     T2[GETMAJIND_YPEN_Ym6]);
 	    }END_FORY1 

	    c2d->transposeY2X_MajorIndex(T2, T);
	    c2d->deallocXYZ(T2);
	}

        FOR_Y1_XPEN{
	    if(neumannLocalY){
	        T[ip] = calcNeumann(T[GETYPIND_XPEN_Ym1],
                                    T[GETYPIND_XPEN_Ym2],
                                    T[GETYPIND_XPEN_Ym3],
                                    T[GETYPIND_XPEN_Ym4],
                                    T[GETYPIND_XPEN_Ym5],
                                    T[GETYPIND_XPEN_Ym6]);
	    } 

            U[ip]  = Y1WallU;
            V[ip]  = 0.0;
            W[ip]  = Y1WallW;
            rhoUP[ip] = rhoP[ip]*Y1WallU;
            rhoVP[ip] = 0.0;
            rhoWP[ip] = rhoP[ip]*Y1WallW;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
        }END_FORY1
    }

    if(bc->bcZ0 == BC::MOVING_ADIABATIC_WALL){

	if(!neumannLocalZ){
	    double *T2, *T3;
	    c2d->allocY(T2);
	    c2d->allocZ(T3);
	    c2d->transposeX2Y_MajorIndex(T, T2);
	    c2d->transposeY2Z_MajorIndex(T2, T3);

	    FOR_Z0_ZPEN_MAJ{
                T3[ip] = calcNeumann(T3[GETMAJIND_ZPEN_Zp1],
                                     T3[GETMAJIND_ZPEN_Zp2],
                                     T3[GETMAJIND_ZPEN_Zp3],
                                     T3[GETMAJIND_ZPEN_Zp4],
                                     T3[GETMAJIND_ZPEN_Zp5],
                                     T3[GETMAJIND_ZPEN_Zp6]);
	    }END_FORZ0

	    c2d->transposeZ2Y_MajorIndex(T3, T2);
	    c2d->transposeY2X_MajorIndex(T2, T);

	    c2d->deallocXYZ(T2);
	    c2d->deallocXYZ(T3);
	}


        FOR_Z0_XPEN{
	    if(neumannLocalZ){
	        T[ip] = calcNeumann(T[GETZPIND_XPEN_Zp1],
                                    T[GETZPIND_XPEN_Zp2],
                                    T[GETZPIND_XPEN_Zp3],
                                    T[GETZPIND_XPEN_Zp4],
                                    T[GETZPIND_XPEN_Zp5],
                                    T[GETZPIND_XPEN_Zp6]);
	    } 

            U[ip]  = Z0WallU;
            V[ip]  = Z0WallV;
            W[ip]  = 0.0;
            rhoUP[ip] = rhoP[ip]*Z0WallU;
            rhoVP[ip] = rhoP[ip]*Z0WallV;
            rhoWP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
        }END_FORZ0
    }

    if(bc->bcZ1 == BC::MOVING_ADIABATIC_WALL){

	if(!neumannLocalZ){
	    double *T2, *T3;
	    c2d->allocY(T2);
	    c2d->allocZ(T3);
	    c2d->transposeX2Y_MajorIndex(T, T2);
	    c2d->transposeY2Z_MajorIndex(T2, T3);

	    FOR_Z1_ZPEN_MAJ{
                T3[ip] = calcNeumann(T3[GETMAJIND_ZPEN_Zm1],
                                     T3[GETMAJIND_ZPEN_Zm2],
                                     T3[GETMAJIND_ZPEN_Zm3],
                                     T3[GETMAJIND_ZPEN_Zm4],
                                     T3[GETMAJIND_ZPEN_Zm5],
                                     T3[GETMAJIND_ZPEN_Zm6]);
 	    }END_FORZ1

	    c2d->transposeZ2Y_MajorIndex(T3, T2);
	    c2d->transposeY2X_MajorIndex(T2, T);

	    c2d->deallocXYZ(T2);
	    c2d->deallocXYZ(T3);
	}

        FOR_Z1_XPEN{
	    if(neumannLocalZ){
	        T[ip] = calcNeumann(T[GETZPIND_XPEN_Zm1],
                                    T[GETZPIND_XPEN_Zm2],
                                    T[GETZPIND_XPEN_Zm3],
                                    T[GETZPIND_XPEN_Zm4],
                                    T[GETZPIND_XPEN_Zm5],
                                    T[GETZPIND_XPEN_Zm6]);
	    } 

            U[ip]  = Z1WallU;
            V[ip]  = Z1WallV;
            W[ip]  = 0.0;
            rhoUP[ip] = rhoP[ip]*Z1WallU;
            rhoVP[ip] = rhoP[ip]*Z1WallV;
            rhoWP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
        }END_FORZ1
    }


    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > preBCHandl Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }



}

void UniformCSolverConservative::postStepBCHandling(){

    if(useTiming) ft1 = MPI_Wtime();

    ////////////////////////////////
    //ADIABATIC AND MOVING WALL BC// 
    ////////////////////////////////

    if(bc->bcX0 == BC::ADIABATIC_WALL || bc->bcX0 == BC::MOVING_ADIABATIC_WALL){
	FOR_X0_XPEN_MAJ{
	    rhok2[ip]  = -ts->dt*contEulerX[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORX0
    }

    if(bc->bcX1 == BC::ADIABATIC_WALL  || bc->bcX1 == BC::MOVING_ADIABATIC_WALL){
	FOR_X1_XPEN_MAJ{
	    rhok2[ip]  = -ts->dt*contEulerX[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORX1
    }   

    if(bc->bcY0 == BC::ADIABATIC_WALL || bc->bcY0 == BC::MOVING_ADIABATIC_WALL){
	FOR_Y0_XPEN_MAJ{
	    rhok2[ip]  = -ts->dt*contEulerY[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORY0
    }

    if(bc->bcY1 == BC::ADIABATIC_WALL || bc->bcY1 == BC::MOVING_ADIABATIC_WALL){
	FOR_Y1_XPEN_MAJ{
	    rhok2[ip]  = -ts->dt*contEulerY[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORY1
    }

    if(bc->bcZ0 == BC::ADIABATIC_WALL || bc->bcZ0 == BC::MOVING_ADIABATIC_WALL){
	FOR_Z0_XPEN_MAJ{
	    rhok2[ip]  = -ts->dt*contEulerZ[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORZ0
    }

    if(bc->bcZ1 == BC::ADIABATIC_WALL || bc->bcZ1 == BC::MOVING_ADIABATIC_WALL){
	FOR_Z1_XPEN_MAJ{
	    rhok2[ip]  = -ts->dt*contEulerZ[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORZ1
    }


    /////////////
    //SPONGE BC//
    /////////////

    if(bc->bcX0 == BC::SPONGE){
	FOR_X0_XPEN_MAJ{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORX0
    }

    if(bc->bcX1 == BC::SPONGE){
	FOR_X1_XPEN_MAJ{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORX1
    }   

    if(bc->bcY0 == BC::SPONGE){
	FOR_Y0_XPEN_MAJ{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORY0
    }

    if(bc->bcY1 == BC::SPONGE){
	FOR_Y1_XPEN_MAJ{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORY1
    }

    if(bc->bcZ0 == BC::SPONGE){
	FOR_Z0_XPEN_MAJ{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORZ0
    }

    if(bc->bcZ1 == BC::SPONGE){
	FOR_Z1_XPEN_MAJ{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORZ1
    }

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > postBCHand Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }

}

void UniformCSolverConservative::updateSponge(){

    if(useTiming) ft1 = MPI_Wtime();

    if(spongeFlag){
	double eps = 1.0/(spg->avgT/ts->dt + 1.0);
	FOR_XYZ_XPEN{
	    spg->spongeRhoAvg[ip]  += eps*(rho1[ip]  - spg->spongeRhoAvg[ip]);	
	    spg->spongeRhoUAvg[ip] += eps*(rhoU1[ip] - spg->spongeRhoUAvg[ip]);	
	    spg->spongeRhoVAvg[ip] += eps*(rhoV1[ip] - spg->spongeRhoVAvg[ip]);	
	    spg->spongeRhoWAvg[ip] += eps*(rhoW1[ip] - spg->spongeRhoWAvg[ip]);	
	    spg->spongeRhoEAvg[ip] += eps*(rhoE1[ip] - spg->spongeRhoEAvg[ip]);	
	    spg->spongeRhoEAvg[ip] = spg->epsP*spg->spongeRhoEAvg[ip] + (1.0 -  spg->epsP)*(spg->spongeP/(ig->gamma-1.0) \
					 + 0.5*(spg->spongeRhoUAvg[ip]*spg->spongeRhoUAvg[ip] + spg->spongeRhoVAvg[ip]*spg->spongeRhoVAvg[ip] \
					 + spg->spongeRhoWAvg[ip]*spg->spongeRhoWAvg[ip])/spg->spongeRhoAvg[ip]);
	}
	
        if(bc->bcX0 == BC::SPONGE){
	    FOR_X0_XPEN{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORX0
        }

        if(bc->bcX1 == BC::SPONGE){
	    FOR_X1_XPEN{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORX1
        }   

        if(bc->bcY0 == BC::SPONGE){
	    FOR_Y0_XPEN{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORY0
        }

        if(bc->bcY1 == BC::SPONGE){
	    FOR_Y1_XPEN{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORY1
        }

        if(bc->bcZ0 == BC::SPONGE){
	    FOR_Z0_XPEN{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORZ0

        }

        if(bc->bcZ1 == BC::SPONGE){
	    FOR_Z1_XPEN{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORZ1
        }

    }

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > updateSpge Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }

};


