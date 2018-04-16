#include "UniformCSolver.hpp"

void UniformCSolver::initializeSolverData(){

    IF_RANK0{
        cout << endl;
        cout << " > Allocating Solver Arrays..." << endl;
        double workSize = 0;
        workSize = 157.0 * (double)N * 8.0;
        cout << " > Need " << workSize/1024.0/1024.0/1024.0 << " Gb of memory required to allocate solver arrays " << endl;
    }

    //3
    c2d->allocX(Ux);
    c2d->allocX(Uy);
    c2d->allocX(Uz);

    //6
    c2d->allocX(Vx);
    c2d->allocX(Vy);
    c2d->allocX(Vz);

    //9
    c2d->allocX(Wx);
    c2d->allocX(Wy);
    c2d->allocX(Wz);

    //12
    c2d->allocX(Uxx);
    c2d->allocX(Uyy);
    c2d->allocX(Uzz);

    //15
    c2d->allocX(Vxx);
    c2d->allocX(Vyy);
    c2d->allocX(Vzz);

    //18
    c2d->allocX(Wxx);
    c2d->allocX(Wyy);
    c2d->allocX(Wzz);

    //21
    c2d->allocX(Uxy);
    c2d->allocX(Uxz);
    c2d->allocX(Uyz);

    //24
    c2d->allocX(Vxy);
    c2d->allocX(Vxz);
    c2d->allocX(Vyz);

    //27
    c2d->allocX(Wxy);
    c2d->allocX(Wxz);
    c2d->allocX(Wyz);

    //30
    c2d->allocX(Tx);
    c2d->allocX(Ty);
    c2d->allocX(Tz);

    //33
    c2d->allocX(Txx);
    c2d->allocX(Tyy);
    c2d->allocX(Tzz);

    //36
    c2d->allocX(contEulerX);
    c2d->allocX(contEulerY);
    c2d->allocX(contEulerZ);

    //39
    c2d->allocX(momXEulerX);
    c2d->allocX(momXEulerY);
    c2d->allocX(momXEulerZ);

    //42
    c2d->allocX(momYEulerX);
    c2d->allocX(momYEulerY);
    c2d->allocX(momYEulerZ);

    //45
    c2d->allocX(momZEulerX);
    c2d->allocX(momZEulerY);
    c2d->allocX(momZEulerZ);

    //48
    c2d->allocX(engyEulerX);
    c2d->allocX(engyEulerY);
    c2d->allocX(engyEulerZ);

    //52
    c2d->allocX(rho1);
    c2d->allocX(rhok);
    c2d->allocX(rhok2);
    c2d->allocX(rho2);

    //56
    c2d->allocX(rhoU1);
    c2d->allocX(rhoUk);
    c2d->allocX(rhoUk2);
    c2d->allocX(rhoU2);

    //60
    c2d->allocX(rhoV1);
    c2d->allocX(rhoVk);
    c2d->allocX(rhoVk2);
    c2d->allocX(rhoV2);

    //64
    c2d->allocX(rhoW1);
    c2d->allocX(rhoWk);
    c2d->allocX(rhoWk2);
    c2d->allocX(rhoW2);
 
    //68
    c2d->allocX(rhoE1);
    c2d->allocX(rhoEk);
    c2d->allocX(rhoEk2);
    c2d->allocX(rhoE2);

    //73 these will be cleared though...
    c2d->allocX(rho0);
    c2d->allocX(U0);
    c2d->allocX(V0);
    c2d->allocX(W0);
    c2d->allocX(p0);

    //76
    c2d->allocX(U);
    c2d->allocY(U2);
    c2d->allocZ(U3);

    //79
    c2d->allocX(V);
    c2d->allocY(V2);
    c2d->allocZ(V3);

    //82
    c2d->allocX(W);
    c2d->allocY(W2);
    c2d->allocZ(W3);

    //85
    c2d->allocX(T);
    c2d->allocY(T2);
    c2d->allocZ(T3);

    //88
    c2d->allocX(p);
    c2d->allocY(p2);
    c2d->allocZ(p3);

    //91
    c2d->allocX(mu);
    c2d->allocX(Amu);
    c2d->allocX(sos);

    //119
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
    c2d->allocY(tempY13);
    c2d->allocY(tempY14);
    c2d->allocY(tempY15);
    c2d->allocY(tempY16);
    c2d->allocY(tempY17);
    c2d->allocY(tempY18);
    c2d->allocY(tempY19);
    c2d->allocY(tempY20);
    c2d->allocY(tempY21);
    c2d->allocY(tempY22);
    c2d->allocY(tempY23);
    c2d->allocY(tempY24);
    c2d->allocY(tempY25);
    c2d->allocY(tempY26);
    c2d->allocY(tempY27);
    c2d->allocY(tempY28);

    //153
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
    c2d->allocZ(tempZ11);
    c2d->allocZ(tempZ12);
    c2d->allocZ(tempZ13);
    c2d->allocZ(tempZ14);
    c2d->allocZ(tempZ15);
    c2d->allocZ(tempZ16);
    c2d->allocZ(tempZ17);
    c2d->allocZ(tempZ18);
    c2d->allocZ(tempZ19);
    c2d->allocZ(tempZ20);
    c2d->allocZ(tempZ21);
    c2d->allocZ(tempZ22);
    c2d->allocZ(tempZ23);
    c2d->allocZ(tempZ24);
    c2d->allocZ(tempZ25);
    c2d->allocZ(tempZ26);
    c2d->allocZ(tempZ27);
    c2d->allocZ(tempZ28);
    c2d->allocZ(tempZ29);
    c2d->allocZ(tempZ30);
    c2d->allocZ(tempZ31);
    c2d->allocZ(tempZ32);
    c2d->allocZ(tempZ33);
    c2d->allocZ(tempZ34);

    //157
    c2d->allocX(tempX1);
    c2d->allocX(tempX2);
    c2d->allocX(tempX3);
    c2d->allocX(tempX4);
}


void UniformCSolver::setInitialConditions(){

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
        }END_FORY0

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

    if(bc->bcY1 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Y1_XPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
       }END_FORY1

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

    if(bc->bcZ0 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Z0_XPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
        }END_FORZ0

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

    if(bc->bcZ1 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Z1_XPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
       }END_FORZ1

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
        }END_FORY0

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

    if(bc->bcY1 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Y1_XPEN_MAJ{
            U[ip]  = Y1WallU;
            V[ip]  = 0.0;
            W[ip]  = Y1WallW;
            rhoU1[ip] = rho1[ip]*Y1WallU;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = rho1[ip]*Y1WallW;
        }END_FORY1


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

    if(bc->bcZ0 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Z0_XPEN_MAJ{
            U[ip]  = Z0WallU;
            V[ip]  = Z0WallV;
            W[ip]  = 0.0;
            rhoU1[ip] = rho1[ip]*Z0WallU;
            rhoV1[ip] = rho1[ip]*Z0WallV;
            rhoW1[ip] = 0.0;
        }END_FORZ0

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

    if(bc->bcZ1 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Z1_XPEN_MAJ{
            U[ip]  = Z1WallU;
            V[ip]  = Z1WallV;
            W[ip]  = 0.0;
            rhoU1[ip] = rho1[ip]*Z1WallU;
            rhoV1[ip] = rho1[ip]*Z1WallV;
            rhoW1[ip] = 0.0;
        }END_FORZ1

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

}

void UniformCSolver::calcDtFromCFL(){
    

    //Calculate the wave speed over the local spacings...
    double *UChar_dx;
    c2d->allocX(UChar_dx);

    FOR_XYZ_XPEN{
	UChar_dx[ip] = (fabs(U[ip]) + sos[ip])/dom->dx + (fabs(V[ip])+sos[ip])/dom->dy + (fabs(W[ip]) + sos[ip])/dom->dz;
    }

    //Get the largest value in the domain
    double lmax_UChar_dx = -100000.0;
    FOR_XYZ_XPEN{
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


}


void UniformCSolver::preStepBCHandling(){

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

	FOR_Y0_XPEN_MAJ{
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


	FOR_Y1_XPEN_MAJ{
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


	FOR_Z0_XPEN{
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


	FOR_Z1_XPEN{
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

        FOR_Y0_XPEN{
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


        FOR_Y1_XPEN{
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

        FOR_Z0_XPEN{
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

        FOR_Z1_XPEN{
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

}


void UniformCSolver::preStepDerivatives(){

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

    ///////////////////
    // X-DERIVATIVES //
    ///////////////////

    //First we'll do all of the X-Direction derivatives since we're in XYZ order

    //Calculate the Euler Components of the equations... 
    FOR_XYZ_XPEN{
	tempX1[ip] = rhoUP[ip]*U[ip] + p[ip];
    	tempX2[ip] = rhoVP[ip]*U[ip];
    	tempX3[ip] = rhoWP[ip]*U[ip];
    	tempX4[ip] = rhoEP[ip]*U[ip] + U[ip]*p[ip];
    }

    derivX->calc1stDerivField(U, Ux); //dU/dx
    derivX->calc1stDerivField(V, Vx); //dV/dx
    derivX->calc1stDerivField(W, Wx); //dW/dx
    derivX->calc1stDerivField(T, Tx); //dT/dx

    derivX->calc2ndDerivField(U, Uxx); //d2U/dx2
    derivX->calc2ndDerivField(V, Vxx); //d2V/dx2
    derivX->calc2ndDerivField(W, Wxx); //d2W/dx2
    derivX->calc2ndDerivField(T, Txx); //d2T/dx2
        
    derivX->calc1stDerivField(rhoUP,  contEulerX); //d(rhoU)/dx
    derivX->calc1stDerivField(tempX1, momXEulerX); //d(rhoU*U + p)/dx
    derivX->calc1stDerivField(tempX2, momYEulerX); //d(rhoU*V)/dx
    derivX->calc1stDerivField(tempX3, momZEulerX); //d(rhoU*W)/dx
    derivX->calc1stDerivField(tempX4, engyEulerX); //d(rhoE*U + U*p)/dx

    c2d->transposeX2Y_MajorIndex(Ux, tempY1); 
    c2d->transposeX2Y_MajorIndex(Vx, tempY2);
    c2d->transposeX2Y_MajorIndex(Wx, tempY3);

    c2d->transposeX2Y_MajorIndex(rhoP,  tempY4);
    c2d->transposeX2Y_MajorIndex(rhoUP, tempY5);
    c2d->transposeX2Y_MajorIndex(rhoVP, tempY6);
    c2d->transposeX2Y_MajorIndex(rhoWP, tempY7);
    c2d->transposeX2Y_MajorIndex(rhoEP, tempY8);

    ///////////////////
    // Y-DERIVATIVES //
    ///////////////////

    FOR_XYZ_YPEN{
        //Now recalculate properties in the new space
	U2[ip] = tempY5[ip]/tempY4[ip];
    	V2[ip] = tempY6[ip]/tempY4[ip];
    	W2[ip] = tempY7[ip]/tempY4[ip];
    	p2[ip] = ig->solvep(tempY4[ip], tempY8[ip], U2[ip], V2[ip], W2[ip]);
    	T2[ip] = ig->solveT(tempY4[ip], p2[ip]);

        //Calculate the stuff for Euler derivatives in new space
    	tempY9[ip]  = tempY5[ip]*V2[ip];
    	tempY10[ip] = tempY6[ip]*V2[ip] + p2[ip];
    	tempY11[ip] = tempY7[ip]*V2[ip];
    	tempY12[ip] = tempY8[ip]*V2[ip] + V2[ip]*p2[ip];
    }


    derivY->calc1stDerivField(U2, tempY13); //dU/dy
    derivY->calc1stDerivField(V2, tempY14); //dV/dy
    derivY->calc1stDerivField(W2, tempY15); //dW/dy

    derivY->calc1stDerivField(tempY1, tempY16); //d2U/dxdy
    derivY->calc1stDerivField(tempY2, tempY17); //d2V/dxdy
    derivY->calc1stDerivField(tempY3, tempY18); //d2W/dxdy
    
    derivY->calc2ndDerivField(U2, tempY19); //d2U/dy2
    derivY->calc2ndDerivField(V2, tempY20); //d2V/dy2
    derivY->calc2ndDerivField(W2, tempY21); //d2W/dy2

    derivY->calc1stDerivField(T2, tempY22); //dT/dy
    derivY->calc2ndDerivField(T2, tempY23); //d2T/dy2

   
    derivY->calc1stDerivField(tempY6, tempY24); //d(rhoV)/dy

    derivY->calc1stDerivField(tempY9, tempY25); //d(rhoU*V)/dy
    derivY->calc1stDerivField(tempY10, tempY26); //d(rhoV*V + p)/dy
    derivY->calc1stDerivField(tempY11, tempY27); //d(rhoW*V)/dy
    derivY->calc1stDerivField(tempY12, tempY28); //d(rhoE*V + V*p)/dy

    //Move things we need to back to X
    c2d->transposeY2X_MajorIndex(tempY13, Uy); 
    c2d->transposeY2X_MajorIndex(tempY14, Vy); 
    c2d->transposeY2X_MajorIndex(tempY15, Wy); 

    c2d->transposeY2X_MajorIndex(tempY16, Uxy);
    c2d->transposeY2X_MajorIndex(tempY17, Vxy);
    c2d->transposeY2X_MajorIndex(tempY18, Wxy);

    c2d->transposeY2X_MajorIndex(tempY19, Uyy);
    c2d->transposeY2X_MajorIndex(tempY20, Vyy);
    c2d->transposeY2X_MajorIndex(tempY21, Wyy);

    c2d->transposeY2X_MajorIndex(tempY22, Ty);
    c2d->transposeY2X_MajorIndex(tempY23, Tyy);

    c2d->transposeY2X_MajorIndex(tempY24, contEulerY);
    c2d->transposeY2X_MajorIndex(tempY25, momXEulerY);
    c2d->transposeY2X_MajorIndex(tempY26, momYEulerY);
    c2d->transposeY2X_MajorIndex(tempY27, momZEulerY);
    c2d->transposeY2X_MajorIndex(tempY28, engyEulerY);

    //Move things we need to move over to Z
    c2d->transposeY2Z_MajorIndex(tempY13, tempZ1); //dU/dy 
    c2d->transposeY2Z_MajorIndex(tempY14, tempZ2); //dV/dy
    c2d->transposeY2Z_MajorIndex(tempY15, tempZ3); //dW/dy

    c2d->transposeY2Z_MajorIndex(tempY1, tempZ4); //dU/dx
    c2d->transposeY2Z_MajorIndex(tempY2, tempZ5); //dV/dx
    c2d->transposeY2Z_MajorIndex(tempY3, tempZ6); //dW/dx

    c2d->transposeY2Z_MajorIndex(tempY4, tempZ7); //rho
    c2d->transposeY2Z_MajorIndex(tempY5, tempZ8); //rhoU
    c2d->transposeY2Z_MajorIndex(tempY6, tempZ9); //rhoV
    c2d->transposeY2Z_MajorIndex(tempY7, tempZ10); //rhoW
    c2d->transposeY2Z_MajorIndex(tempY8, tempZ11); //rhoE

    ///////////////////
    // Z-DERIVATIVES //
    ///////////////////

    //Now recalculate properties in the new space
    FOR_XYZ_ZPEN{
        U3[ip] =  tempZ8[ip]/tempZ7[ip];
        V3[ip] =  tempZ9[ip]/tempZ7[ip];
        W3[ip] = tempZ10[ip]/tempZ7[ip];
        p3[ip] = ig->solvep(tempZ7[ip], tempZ11[ip], U3[ip], V3[ip], W3[ip]);
        T3[ip] = ig->solveT(tempZ7[ip], p3[ip]);

        //Calculate the stuff for the Euler Derivatives
        tempZ12[ip] =  tempZ8[ip]*W3[ip];
        tempZ13[ip] =  tempZ9[ip]*W3[ip];
        tempZ14[ip] = tempZ10[ip]*W3[ip] + p3[ip];
        tempZ15[ip] = tempZ11[ip]*W3[ip] + W3[ip]*p3[ip];
    }

//TODO LEFT OFF MOVING TO MPI HERE!!!
   
    //Calculate the Z direction derivatives
    derivZ->calc1stDerivField(U3, tempZ16); //dU/dz
    derivZ->calc1stDerivField(V3, tempZ17); //dV/dz
    derivZ->calc1stDerivField(W3, tempZ18); //dW/dz
    derivZ->calc1stDerivField(T3, tempZ19); //dT/dz

    derivZ->calc2ndDerivField(U3, tempZ20); //d2U/dz2
    derivZ->calc2ndDerivField(V3, tempZ21); //d2V/dz2
    derivZ->calc2ndDerivField(W3, tempZ22); //d2W/dz2
    derivZ->calc2ndDerivField(T3, tempZ23); //d2T/dz2

    derivZ->calc1stDerivField(tempZ1, tempZ24); //d2U/dydz
    derivZ->calc1stDerivField(tempZ2, tempZ25); //d2V/dydz
    derivZ->calc1stDerivField(tempZ3, tempZ26); //d2W/dydz

    derivZ->calc1stDerivField(tempZ4, tempZ27); //d2U/dxdz 
    derivZ->calc1stDerivField(tempZ5, tempZ28); //d2V/dxdz 
    derivZ->calc1stDerivField(tempZ6, tempZ29); //d2W/dxdz 

    derivZ->calc1stDerivField(tempZ12, tempZ30); //d(rhoU*W)/dz
    derivZ->calc1stDerivField(tempZ13, tempZ31); //d(rhoV*W)/dz
    derivZ->calc1stDerivField(tempZ14, tempZ32); //d(rhoW*W + p)/dz
    derivZ->calc1stDerivField(tempZ15, tempZ33); //d(rhoE*W + W*p)/dz
    derivZ->calc1stDerivField(tempZ10, tempZ34); //d(rhoW)/dz
    
    //Move the things we need back to X
    c2d->transposeZ2Y_MajorIndex(tempZ16, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, Uz);

    c2d->transposeZ2Y_MajorIndex(tempZ17, tempY1);
    c2d->transposeY2X_MajorIndex(tempY1, Vz);
 
    c2d->transposeZ2Y_MajorIndex(tempZ18, tempY1);
    c2d->transposeY2X_MajorIndex(tempY1, Wz);
 
    c2d->transposeZ2Y_MajorIndex(tempZ19, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, Ty);

    c2d->transposeZ2Y_MajorIndex(tempZ20, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, Uzz);

    c2d->transposeZ2Y_MajorIndex(tempZ21, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, Vzz);

    c2d->transposeZ2Y_MajorIndex(tempZ22, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, Wzz);

    c2d->transposeZ2Y_MajorIndex(tempZ23, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, Tzz);

    c2d->transposeZ2Y_MajorIndex(tempZ24, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, Uyz);

    c2d->transposeZ2Y_MajorIndex(tempZ25, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, Vyz);

    c2d->transposeZ2Y_MajorIndex(tempZ26, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, Wyz);

    c2d->transposeZ2Y_MajorIndex(tempZ27, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, Uxz);

    c2d->transposeZ2Y_MajorIndex(tempZ28, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, Vxz);

    c2d->transposeZ2Y_MajorIndex(tempZ29, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, Wxz);

    c2d->transposeZ2Y_MajorIndex(tempZ34, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, contEulerZ);

    c2d->transposeZ2Y_MajorIndex(tempZ30, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, momXEulerZ);

    c2d->transposeZ2Y_MajorIndex(tempZ31, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, momYEulerZ);

    c2d->transposeZ2Y_MajorIndex(tempZ32, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, momZEulerZ);

    c2d->transposeZ2Y_MajorIndex(tempZ33, tempY1); 
    c2d->transposeY2X_MajorIndex(tempY1, engyEulerZ);

}


void UniformCSolver::solveContinuity(){

    double *rhoP;
    if(rkStep == 1){
        rhoP = rho1;
    }else{
        rhoP = rhok;
    }
	
    double spgSource; 

    FOR_XYZ_XPEN{

	if(spongeFlag)
	    spgSource = calcSpongeSource(rhoP[ip], spg->spongeRhoAvg[ip], spg->sigma[ip]);
	else
	    spgSource = 0.0;
		
	rhok2[ip]  = ts->dt*(-contEulerX[ip] - contEulerY[ip] - contEulerZ[ip] + spgSource);
    }

}

void UniformCSolver::solveXMomentum(){


    double *rhoUP;
    if(rkStep == 1){
        rhoUP = rhoU1;
    }else{
        rhoUP = rhoUk;
    }

    double MuX, MuY, MuZ, spgSource;
    FOR_XYZ_XPEN{

	if(spongeFlag)
	    spgSource = calcSpongeSource(rhoUP[ip], spg->spongeRhoUAvg[ip], spg->sigma[ip]);
	else
	    spgSource = 0.0;

	//Calculate the viscosity derivatives
	MuX = Amu[ip]*Tx[ip];
	MuY = Amu[ip]*Ty[ip];
	MuZ = Amu[ip]*Tz[ip];

	//Viscous Terms
        rhoUk2[ip]  = mu[ip]*((4.0/3.0)*Uxx[ip] + Uyy[ip] + Uzz[ip] + (1.0/3.0)*Vxy[ip] + (1.0/3.0)*Wxz[ip]);
	rhoUk2[ip] += (4.0/3.0)*MuX*(Ux[ip] - 0.5*Vy[ip] - 0.5*Wz[ip]) + MuY*(Uy[ip] + Vx[ip]) + MuZ*(Wx[ip] + Uz[ip]);

	//Euler Terms
	rhoUk2[ip] += -momXEulerX[ip] -momXEulerY[ip] -momXEulerZ[ip] + spgSource;
	rhoUk2[ip] *= ts->dt;
    }


}

void UniformCSolver::solveYMomentum(){

    double *rhoVP;
    if(rkStep == 1){
        rhoVP = rhoV1;
    }else{
        rhoVP = rhoVk;
    }

    double MuY, MuX, MuZ, spgSource;

    FOR_XYZ_XPEN{ 

        if(spongeFlag)
            spgSource = calcSpongeSource(rhoVP[ip], spg->spongeRhoVAvg[ip], spg->sigma[ip]);
	else
	    spgSource = 0.0;

		
	MuX = Amu[ip]*Tx[ip];
	MuY = Amu[ip]*Ty[ip];
  	MuZ = Amu[ip]*Tz[ip];

        //Viccous Terms 
	rhoVk2[ip]  = mu[ip]*((4.0/3.0)*Vyy[ip] + Vxx[ip] + Vzz[ip] + (1.0/3.0)*Uxy[ip] + (1.0/3.0)*Wyz[ip]);
	rhoVk2[ip] += (4.0/3.0)*MuY*(Vy[ip] - 0.5*Ux[ip] - 0.5*Wz[ip]) + MuX*(Uy[ip] + Vx[ip]) + MuZ*(Wy[ip] + Vz[ip]);
	//Euler Terms
	rhoVk2[ip] += -momYEulerX[ip] -momYEulerY[ip] -momYEulerZ[ip] + spgSource;

        rhoVk2[ip] *= ts->dt;
   }

}

void UniformCSolver::solveZMomentum(){

    double *rhoWP;
    if(rkStep == 1){
        rhoWP = rhoW1;
    }else{
        rhoWP = rhoWk;
    }

    double MuY, MuX, MuZ, spgSource;


    FOR_XYZ_XPEN{

        if(spongeFlag)
            spgSource = calcSpongeSource(rhoWP[ip], spg->spongeRhoWAvg[ip], spg->sigma[ip]);
        else
	    spgSource = 0.0;

	MuZ = Amu[ip]*Tz[ip];
	MuX = Amu[ip]*Tx[ip];
	MuY = Amu[ip]*Ty[ip];


    	//Viscous Terms
        rhoWk2[ip]  = mu[ip]*((4.0/3.0)*Wzz[ip] + Wyy[ip] + Wxx[ip] + (1.0/3.0)*Uxz[ip] + (1.0/3.0)*Vyz[ip]);
	rhoWk2[ip] += (4.0/3.0)*MuZ*(Wz[ip] - 0.5*Ux[ip] - 0.5*Vy[ip]) + MuX*(Wx[ip] + Uz[ip]) + MuY*(Wy[ip] + Vz[ip]);

	//Euler Terms
	rhoWk2[ip] += -momZEulerX[ip] -momZEulerY[ip] -momZEulerZ[ip] + spgSource;

        rhoWk2[ip] *= ts->dt;
    }

}


void UniformCSolver::solveEnergy(){


    double *rhoEP;
    if(rkStep == 1){
        rhoEP = rhoE1;
    }else{
        rhoEP = rhoEk;
    }

    double qtemp, vtemp1, vtemp2, engyEuler;
    double MuX, MuY, MuZ, spgSource;

    FOR_XYZ_XPEN{

        if(spongeFlag)
            spgSource = calcSpongeSource(rhoEP[ip], spg->spongeRhoEAvg[ip], spg->sigma[ip]);
        else
            spgSource = 0.0;

        MuX = Amu[ip]*Tx[ip];
	MuY = Amu[ip]*Ty[ip];
	MuZ = Amu[ip]*Tz[ip];

    	//Heat Transfer Terms
	qtemp   =  ig->cp/ig->Pr*(MuX*Tx[ip]     + MuY*Ty[ip]     +  MuZ*Tz[ip] + 
		     	 	  mu[ip]*Txx[ip] + mu[ip]*Tyy[ip] + mu[ip]*Tzz[ip]);



    	//Viscous Energy terms w/o viscosity derivatives...
	vtemp1  = mu[ip]*(U[ip]*((4.0/3.0)*Uxx[ip] + Uyy[ip] + Uzz[ip]) + 
		          V[ip]*(Vxx[ip] + (4.0/3.0)*Vyy[ip] + Vzz[ip]) + 
		          W[ip]*(Wxx[ip] + Wyy[ip] + (4.0/3.0)*Wzz[ip]) + 
		(4.0/3.0)*(Ux[ip]*Ux[ip] + Vy[ip]*Vy[ip] + Wz[ip]*Wz[ip]) +
		           Uy[ip]*Uy[ip] + Uz[ip]*Uz[ip] + 
		           Vx[ip]*Vx[ip] + Vz[ip]*Vz[ip] + 
		           Wx[ip]*Wx[ip] + Wy[ip]*Wy[ip] -
	        (4.0/3.0)*(Ux[ip]*Vy[ip] + Ux[ip]*Wz[ip] + Vy[ip]*Wz[ip]) +
		      2.0*(Uy[ip]*Vx[ip] + Uz[ip]*Wx[ip] + Vz[ip]*Wy[ip]) +
		(1.0/3.0)*(U[ip]*Vxy[ip] + U[ip]*Wxz[ip] + V[ip]*Uxy[ip]) +
		(1.0/3.0)*(V[ip]*Wyz[ip] + W[ip]*Uxz[ip] + W[ip]*Vyz[ip]));


    	//Viscous Energy terms w/ viscosity derivatives...
	vtemp2  =  (4.0/3.0)*(U[ip]*MuX*Ux[ip] + V[ip]*MuY*Vy[ip] + W[ip]*MuZ*Wz[ip]) -
		   (2.0/3.0)* U[ip]*MuX*(Vy[ip] + Wz[ip]) -
		   (2.0/3.0)* V[ip]*MuY*(Ux[ip] + Wz[ip]) -
		   (2.0/3.0)* W[ip]*MuZ*(Ux[ip] + Vy[ip]) +
			      U[ip]*MuY*(Uy[ip] + Vx[ip]) +
			      U[ip]*MuZ*(Uz[ip] + Wx[ip]) + 
			      V[ip]*MuX*(Uy[ip] + Vx[ip]) +
			      V[ip]*MuZ*(Vz[ip] + Wy[ip]) +
			      W[ip]*MuX*(Uz[ip] + Wx[ip]) +
			      W[ip]*MuY*(Vz[ip] + Wy[ip]);


    	//Euler terms
	engyEuler  = -engyEulerX[ip] - engyEulerY[ip] - engyEulerZ[ip] + spgSource;


	//Put it all together...
	rhoEk2[ip] = ts->dt*(qtemp + vtemp1 + vtemp2 + engyEuler + spgSource);
    }

}

void UniformCSolver::postStepBCHandling(){

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

    ////////////////////////////////
    //ADIABATIC AND MOVING WALL BC// 
    ////////////////////////////////

    if(bc->bcX0 == BC::ADIABATIC_WALL || bc->bcX0 == BC::MOVING_ADIABATIC_WALL){
	FOR_X0_XPEN_MAJ{
	    rhok2[ip]  = -ts->dt*contEulerX[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    //rhoEk2[ip] = -ts->dt*(engyEulerX[ip] - (mu[ip]/ig->Pr/(ig->gamma-1.0))*Txx[ip]);
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORX0
    }

    if(bc->bcX1 == BC::ADIABATIC_WALL  || bc->bcX1 == BC::MOVING_ADIABATIC_WALL){
	FOR_X1_XPEN_MAJ{
	    rhok2[ip]  = -ts->dt*contEulerX[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    //rhoEk2[ip] = -ts->dt*(engyEulerX[ip] - (mu[ip]/ig->Pr/(ig->gamma-1.0))*Txx[ip]);
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORX1
    }   

    if(bc->bcY0 == BC::ADIABATIC_WALL || bc->bcY0 == BC::MOVING_ADIABATIC_WALL){
	FOR_Y0_XPEN_MAJ{
	    rhok2[ip]  = -ts->dt*contEulerY[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    //rhoEk2[ip] = -ts->dt*(engyEulerY[ip] - (mu[ip]/ig->Pr/(ig->gamma-1.0))*Tyy[ip]);
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORY0
    }

    if(bc->bcY1 == BC::ADIABATIC_WALL || bc->bcY1 == BC::MOVING_ADIABATIC_WALL){
	FOR_Y1_XPEN_MAJ{
	    rhok2[ip]  = -ts->dt*contEulerY[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    //rhoEk2[ip] = -ts->dt*(engyEulerY[ip] - (mu[ip]/ig->Pr/(ig->gamma-1.0))*Tyy[ip]);
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORY1
    }

    if(bc->bcZ0 == BC::ADIABATIC_WALL || bc->bcZ0 == BC::MOVING_ADIABATIC_WALL){
	FOR_Z0_XPEN_MAJ{
	    rhok2[ip]  = -ts->dt*contEulerZ[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    //rhoEk2[ip] = -ts->dt*(engyEulerZ[ip] - (mu[ip]/ig->Pr/(ig->gamma-1.0))*Tzz[ip]);
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORZ0
    }

    if(bc->bcZ1 == BC::ADIABATIC_WALL || bc->bcZ1 == BC::MOVING_ADIABATIC_WALL){
	FOR_Z1_XPEN_MAJ{
	    rhok2[ip]  = -ts->dt*contEulerZ[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    //rhoEk2[ip] = -ts->dt*(engyEulerZ[ip] - (mu[ip]/ig->Pr/(ig->gamma-1.0))*Tzz[ip]);
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

}


void UniformCSolver::filterConservedData(){

    //Need to do round robin of filtering directions...
    if(timeStep%ts->filterStep == 0){

        //Advance the filtering time step       
        filterTimeStep++;

        //Going to try and be cute to minimize dmemory allocation
        if(filterTimeStep%3 == 1){

            //Here we'll do X->Y->Z     

                    filtX->filterField(rho2,  rho1);
		    filtX->filterField(rhoU2, rhoU1);
		    filtX->filterField(rhoV2, rhoV1);
		    filtX->filterField(rhoW2, rhoW1);
		    filtX->filterField(rhoE2, rhoE1);

		    c2d->transposeX2Y_MajorIndex(rho1,  tempY1);
		    c2d->transposeX2Y_MajorIndex(rhoU1, tempY2);
		    c2d->transposeX2Y_MajorIndex(rhoV1, tempY3);
		    c2d->transposeX2Y_MajorIndex(rhoW1, tempY4);
		    c2d->transposeX2Y_MajorIndex(rhoE1, tempY5);

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

		    c2d->transposeZ2Y_MajorIndex(tempZ6,  tempY1);
		    c2d->transposeZ2Y_MajorIndex(tempZ7,  tempY2);
		    c2d->transposeZ2Y_MajorIndex(tempZ8,  tempY3);
		    c2d->transposeZ2Y_MajorIndex(tempZ9,  tempY4);
		    c2d->transposeZ2Y_MajorIndex(tempZ10, tempY5);

		    c2d->transposeY2X_MajorIndex(tempY1, rho1);
		    c2d->transposeY2X_MajorIndex(tempY2, rhoU1);
		    c2d->transposeY2X_MajorIndex(tempY3, rhoV1);
		    c2d->transposeY2X_MajorIndex(tempY4, rhoW1);
		    c2d->transposeY2X_MajorIndex(tempY5, rhoE1);

        }else if(filterTimeStep%3 == 2){

            //Here we'll do Y->Z->X     

	    c2d->transposeX2Y_MajorIndex(rho2,  tempY1);
	    c2d->transposeX2Y_MajorIndex(rhoU2, tempY2);
	    c2d->transposeX2Y_MajorIndex(rhoV2, tempY3);
	    c2d->transposeX2Y_MajorIndex(rhoW2, tempY4);
	    c2d->transposeX2Y_MajorIndex(rhoE2, tempY5);

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

	    c2d->transposeZ2Y_MajorIndex(tempZ6,  tempY1);
	    c2d->transposeZ2Y_MajorIndex(tempZ7,  tempY2);
	    c2d->transposeZ2Y_MajorIndex(tempZ8,  tempY3);
	    c2d->transposeZ2Y_MajorIndex(tempZ9,  tempY4);
	    c2d->transposeZ2Y_MajorIndex(tempZ10, tempY5);

	    c2d->transposeY2X_MajorIndex(tempY1,  rho2);
	    c2d->transposeY2X_MajorIndex(tempY2,  rhoU2);
	    c2d->transposeY2X_MajorIndex(tempY3,  rhoV2);
	    c2d->transposeY2X_MajorIndex(tempY4,  rhoW2);
	    c2d->transposeY2X_MajorIndex(tempY5,  rhoE2);

	    filtX->filterField(rho2,  rho1);
	    filtX->filterField(rhoU2, rhoU1);
	    filtX->filterField(rhoV2, rhoV1);
	    filtX->filterField(rhoW2, rhoW1);
	    filtX->filterField(rhoE2, rhoE1);

        }else{

	    c2d->transposeX2Y_MajorIndex(rho1,  tempY1);
	    c2d->transposeX2Y_MajorIndex(rhoU1, tempY2);
	    c2d->transposeX2Y_MajorIndex(rhoV1, tempY3);
	    c2d->transposeX2Y_MajorIndex(rhoW1, tempY4);
	    c2d->transposeX2Y_MajorIndex(rhoE1, tempY5);

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

	    c2d->transposeY2X_MajorIndex(tempY1,  rho2);
	    c2d->transposeY2X_MajorIndex(tempY2,  rhoU2);
	    c2d->transposeY2X_MajorIndex(tempY3,  rhoV2);
	    c2d->transposeY2X_MajorIndex(tempY4,  rhoW2);
	    c2d->transposeY2X_MajorIndex(tempY5,  rhoE2);

	    filtX->filterField(rho2,  rho1);
	    filtX->filterField(rhoU2, rhoU1);
	    filtX->filterField(rhoV2, rhoV1);
	    filtX->filterField(rhoW2, rhoW1);
	    filtX->filterField(rhoE2, rhoE1);

	    c2d->transposeX2Y_MajorIndex(rho1,  tempY1);
	    c2d->transposeX2Y_MajorIndex(rhoU1, tempY2);
	    c2d->transposeX2Y_MajorIndex(rhoV1, tempY3);
	    c2d->transposeX2Y_MajorIndex(rhoW1, tempY4);
	    c2d->transposeX2Y_MajorIndex(rhoE1, tempY5);

	    filtY->filterField(tempY1, tempY6);
	    filtY->filterField(tempY2, tempY7);
	    filtY->filterField(tempY3, tempY8);
	    filtY->filterField(tempY4, tempY9);
	    filtY->filterField(tempY5, tempY10);

	    c2d->transposeY2X_MajorIndex(tempY6,  rho1);
	    c2d->transposeY2X_MajorIndex(tempY7,  rhoU1);
	    c2d->transposeY2X_MajorIndex(tempY8,  rhoV1);
	    c2d->transposeY2X_MajorIndex(tempY9,  rhoW1);
	    c2d->transposeY2X_MajorIndex(tempY10, rhoE1);

        }
 
    //If not filtering, need to copy the solution over to the *1 variables
    }else{
	    memcpy(rho1,  rho2,  sizeof(double)*pxSize[0]*pxSize[1]*pxSize[2]);
	    memcpy(rhoU1, rhoU2, sizeof(double)*pxSize[0]*pxSize[1]*pxSize[2]);
	    memcpy(rhoV1, rhoV2, sizeof(double)*pxSize[0]*pxSize[1]*pxSize[2]);
	    memcpy(rhoW1, rhoW2, sizeof(double)*pxSize[0]*pxSize[1]*pxSize[2]);
	    memcpy(rhoE1, rhoE2, sizeof(double)*pxSize[0]*pxSize[1]*pxSize[2]);

    }



};


void UniformCSolver::updateNonConservedData(){

    if(!rkLast){

	FOR_XYZ_XPEN{
	    U[ip]   = ig->solveU(rhok[ip], rhoUk[ip]);
	    V[ip]   = ig->solveU(rhok[ip], rhoVk[ip]);
	    W[ip]   = ig->solveU(rhok[ip], rhoWk[ip]);
	    p[ip]   = ig->solvep(rhok[ip], rhoEk[ip], U[ip], V[ip], W[ip]);
	    T[ip]   = ig->solveT(rhok[ip], p[ip]);
	    mu[ip]  = ig->solveMu(T[ip]);
	    Amu[ip] = ig->solveAmu(T[ip]);
	    sos[ip] = ig->solveSOS(rhok[ip], p[ip]);
	}

    }else if(rkLast){

	FOR_XYZ_XPEN{
	    U[ip]   = ig->solveU(rho1[ip], rhoU1[ip]);
	    V[ip]   = ig->solveU(rho1[ip], rhoV1[ip]);
	    W[ip]   = ig->solveU(rho1[ip], rhoW1[ip]);
	    p[ip]   = ig->solvep(rho1[ip], rhoE1[ip], U[ip], V[ip], W[ip]);
	    T[ip]   = ig->solveT(rho1[ip], p[ip]);
	    mu[ip]  = ig->solveMu(T[ip]);
	    Amu[ip] = ig->solveAmu(T[ip]);
	    sos[ip] = ig->solveSOS(rho1[ip], p[ip]);
	}
    }

}

void UniformCSolver::updateSponge(){

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

};
/*
void UniformCSolver::calcTurbulenceQuantities(){

    if(useTiming) tic();

    #pragma omp parallel for
    FOR_XYZ{
	double Sij[3][3];
	Sij[0][0] = 0.5*(Ux[ip] + Ux[ip]) - (1.0/3.0)*(Ux[ip] + Vy[ip] + Wz[ip]);
	Sij[1][1] = 0.5*(Vy[ip] + Vy[ip]) - (1.0/3.0)*(Ux[ip] + Vy[ip] + Wz[ip]);
	Sij[2][2] = 0.5*(Wz[ip] + Wz[ip]) - (1.0/3.0)*(Ux[ip] + Vy[ip] + Wz[ip]);
	Sij[0][1] = 0.5*(Uy[ip] + Vx[ip]);
	Sij[1][0] = 0.5*(Uy[ip] + Vx[ip]);
	Sij[0][2] = 0.5*(Uz[ip] + Wx[ip]);
	Sij[2][0] = 0.5*(Uz[ip] + Wx[ip]);
	Sij[1][2] = 0.5*(Vz[ip] + Wy[ip]);
	Sij[2][1] = 0.5*(Vz[ip] + Wy[ip]);

	turbdiss[ip]  = Sij[0][0]*Ux[ip];
	turbdiss[ip] += Sij[0][1]*Uy[ip];
	turbdiss[ip] += Sij[0][2]*Uz[ip];
	turbdiss[ip] += Sij[1][0]*Vx[ip];
	turbdiss[ip] += Sij[1][1]*Vy[ip];
	turbdiss[ip] += Sij[1][2]*Vz[ip];
	turbdiss[ip] += Sij[2][0]*Wx[ip];
	turbdiss[ip] += Sij[2][1]*Wy[ip];
	turbdiss[ip] += Sij[2][2]*Wz[ip];
	turbdiss[ip] *= 2*mu[ip];
	
	uprime2[ip]    = (U[ip]*U[ip] + V[ip]*V[ip] + W[ip]*W[ip])/3.0; 
	uvar[ip]       = (U[ip]*U[ip] + V[ip]*V[ip] + W[ip]*W[ip]); 
	uiprime2[ip]   = (Ux[ip]*Ux[ip] + Vy[ip]*Vy[ip] + Wz[ip]*Wz[ip])/3.0; 
	kineticEng[ip] = rho1[ip]*uvar[ip];

    }

    double taylor2_denom = 0;
    double meanRho = 0.0;
    double meanMu = 0.0;
    double meanTurbDiss = 0.0;
    double meanSOS = 0.0; 

    double turbMach = 0.0;
    double uprime = 0.0;
    double uiprime = 0.0;
    double urms = 0.0;
    double kolNu = 0.0;
    double meanKineticEng = 0.0;
 
    double rhoprime = 0.0;
    double dilprime = 0.0;
    double vortprime = 0.0;

    #pragma omp parallel for reduction(+:meanRho, meanMu, meanTurbDiss, uprime, uiprime, urms, meanSOS, taylor2_denom, meanKineticEng) 
    FOR_XYZ{
	meanRho += rho1[ip]/((double)(Nx*Ny*Nz));
	meanMu  += mu[ip]/((double)(Nx*Ny*Nz));
	meanTurbDiss += turbdiss[ip]/((double)(Nx*Ny*Nz));
	uprime += uprime2[ip]/((double)(Nx*Ny*Nz));
	uiprime += uiprime2[ip]/((double)(Nx*Ny*Nz));
	urms   += uvar[ip]/((double)(Nx*Ny*Nz));
	meanSOS += sos[ip]/((double)(Nx*Ny*Nz));
	taylor2_denom += (Ux[ip]*Ux[ip])/((double)(Nx*Ny*Nz));
	meanKineticEng += kineticEng[ip]/((double)(Nx*Ny*Nz))/2.0;
    }

    #pragma omp parallel for reduction(+:rhoprime, dilprime, vortprime)
    FOR_XYZ{
	rhoprime += ((rho1[ip] - meanRho)*(rho1[ip] - meanRho))/((double)(Nx*Ny*Nz));
	dilprime += ((Ux[ip] + Vy[ip] + Wz[ip])*(Ux[ip] + Vy[ip] + Wz[ip]))/((double)(Nx*Ny*Nz));

	double wx, wy, wz;
        wx = Wy[ip] - Vz[ip];
        wy = Uz[ip] - Wx[ip];
        wz = Vx[ip] - Uy[ip];
	vortprime += ((wx*wx + wy*wy + wz*wz))/((double)(Nx*Ny*Nz));
	
    }

    double enstrophy2Mu = meanMu*vortprime;

    rhoprime  = sqrt(rhoprime);
    dilprime  = sqrt(dilprime);
    vortprime = sqrt(vortprime);

    double meanNu = meanMu/meanRho; 
    meanTurbDiss /= meanRho;
    uprime = sqrt(uprime);
    uiprime = sqrt(uiprime);
    urms = sqrt(urms);

    turbMach = urms/meanSOS;

    double taylorMicro2 = uprime/uiprime; //this is the correct one for Samtaney paper
    double taylorMicro  = sqrt(15.0 * meanNu / meanTurbDiss)*urms;
    double taylorReyn  = uprime * taylorMicro / meanNu;
    double taylorReyn2 = urms*urms*sqrt(5.0/(meanNu*meanTurbDiss));
    double taylorReyn3  = uprime * taylorMicro2 / meanNu; //correct one for Samtaney paper
    double k0 = 8.0;
    double tau = sqrt(2.0*M_PI)/k0/uprime;


    //Kolmogorov Scale
    kolNu = pow(pow(meanNu,3.0)/meanTurbDiss,0.25);

    cout << " taylorReyn3: "  << taylorReyn3 << endl;
    cout << " turbMach: "  << turbMach << endl;
    cout << " taylorMicro2: " << taylorMicro2 << endl;
    cout << " meanTurbDiss: " << meanTurbDiss << endl;
    cout << " enstropy2Mu: " << enstrophy2Mu << endl;
    cout << " uprime: "      << uprime << endl;
    cout << " urms: "        << urms << endl;
    cout << " meanRho: "     << meanRho << endl;
    cout << " meanMu: "      << meanMu << endl;
    cout << " KolmogorovNu: " << kolNu << endl;
    cout << " rhoprime: " << rhoprime << endl;
    cout << " dilprime: " << dilprime << endl;
    cout << " vortprime: " << vortprime << endl;
    cout << " tau: " << tau << endl;
    cout << " t/tau: " << time/tau << endl;
    cout << " meanKineticEng: " << meanKineticEng << endl << endl;

        ofstream outfile;
        outfile.precision(17);
        string outputFileName;
        outputFileName = "turbdata.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
        outfile << time << " " << tau << " " << taylorReyn3 << " " << turbMach << " ";
	outfile << meanTurbDiss << " " << uprime << " " << kolNu << " ";
	outfile << rhoprime << " " << dilprime << " " << vortprime << " " << meanKineticEng << " ";
	outfile << endl;
        outfile.close();

    if(useTiming){
        cout << " > calcTrbDat Timing: ";
        toc();
    }



};

void UniformCSolver::calcTaylorGreenQuantities(){


    if(useTiming) tic();

    double enstrophySum = 0.0;
    double enstrophySum2Mu = 0.0;
    double turbDissSum = 0.0;
    double kineticEngSum = 0.0;
    double rhoprime = 0.0;
    double dilprime = 0.0;
    double meanMu   = 0.0;

    #pragma omp parallel for reduction(+:enstrophySum, enstrophySum2Mu, turbDissSum, kineticEngSum, dilprime, meanMu)
    FOR_XYZ{
	double Sij[3][3];
	Sij[0][0] = 0.5*(Ux[ip] + Ux[ip]) - (1.0/3.0)*(Ux[ip] + Vy[ip] + Wz[ip]);
	Sij[1][1] = 0.5*(Vy[ip] + Vy[ip]) - (1.0/3.0)*(Ux[ip] + Vy[ip] + Wz[ip]);
	Sij[2][2] = 0.5*(Wz[ip] + Wz[ip]) - (1.0/3.0)*(Ux[ip] + Vy[ip] + Wz[ip]);
	Sij[0][1] = 0.5*(Uy[ip] + Vx[ip]);
	Sij[1][0] = 0.5*(Uy[ip] + Vx[ip]);
	Sij[0][2] = 0.5*(Uz[ip] + Wx[ip]);
	Sij[2][0] = 0.5*(Uz[ip] + Wx[ip]);
	Sij[1][2] = 0.5*(Vz[ip] + Wy[ip]);
	Sij[2][1] = 0.5*(Vz[ip] + Wy[ip]);

	turbdiss[ip]  = Sij[0][0]*Ux[ip];
	turbdiss[ip] += Sij[0][1]*Uy[ip];
	turbdiss[ip] += Sij[0][2]*Uz[ip];
	turbdiss[ip] += Sij[1][0]*Vx[ip];
	turbdiss[ip] += Sij[1][1]*Vy[ip];
	turbdiss[ip] += Sij[1][2]*Vz[ip];
	turbdiss[ip] += Sij[2][0]*Wx[ip];
	turbdiss[ip] += Sij[2][1]*Wy[ip];
	turbdiss[ip] += Sij[2][2]*Wz[ip];
	turbdiss[ip] *= 2*mu[ip];

        uvar[ip]       = (U[ip]*U[ip] + V[ip]*V[ip] + W[ip]*W[ip]);
	kineticEng[ip] = rho1[ip]*uvar[ip];

	double wx, wy, wz;
        wx = Wy[ip] - Vz[ip];
        wy = Uz[ip] - Wx[ip];
        wz = Vx[ip] - Uy[ip];
	enstrophySum += (1.0/2.0)*rho1[ip]*(wx*wx + wy*wy + wz*wz);
	enstrophySum2Mu += mu[ip]*rho1[ip]*(wx*wx + wy*wy + wz*wz);
	
	turbDissSum += turbdiss[ip];
	kineticEngSum += kineticEng[ip]/2.0;

	dilprime += ((Ux[ip] + Vy[ip] + Wz[ip])*(Ux[ip] + Vy[ip] + Wz[ip]))/((double)(Nx*Ny*Nz));

	meanMu += mu[ip]/((double)(Nx*Ny*Nz));
    }

    dilprime  = sqrt(dilprime);


    cout << " > kineticEngSum: " << kineticEngSum <<  endl;
    cout << " > turbDissSum:   " << turbDissSum <<  endl;
    cout << " > enstrophySum:  " << enstrophySum <<  endl;
    cout << " > enstrophySum*2mu :  " << enstrophySum2Mu <<  endl;
    cout << " > meanMu:        " << meanMu << endl;
    cout << " > dilprime:      " << dilprime << endl << endl;

        ofstream outfile;
        outfile.precision(17);
        string outputFileName;
        outputFileName = "taylorgreen.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
        outfile << time << " " << kineticEngSum << " " << turbDissSum << " " << enstrophySum << " " << enstrophySum2Mu << " " << meanMu << " " << dilprime;
	outfile << endl;
        outfile.close();

    if(useTiming){
        cout << " > calcTayGrn Timing: ";
        toc();
    }



};

void UniformCSolver::shearLayerInfoCalc(){


    if(useTiming) tic();

    if(timeStep%25 == 0){

	double rhoAvg[Ny], uAvg[Ny], vAvg[Ny], wAvg[Ny], pAvg[Ny], tAvg[Ny];
	double rhoprime2[Ny], uprime2[Ny], vprime2[Ny], wprime2[Ny], pprime2[Ny], tprime2[Ny];
	double uvprime[Ny], uwprime[Ny], vwprime[Ny];

	//Calculate the averages...	
	#pragma omp parallel for
	FOR_Y{
	    rhoAvg[j] = 0.0;
	    uAvg[j] = 0.0;
	    vAvg[j] = 0.0;
	    wAvg[j] = 0.0;
	    pAvg[j] = 0.0;
	    tAvg[j] = 0.0;

	    rhoprime2[j] = 0.0;
	    uprime2[j] = 0.0;
	    vprime2[j] = 0.0;
	    wprime2[j] = 0.0;
	    pprime2[j] = 0.0;
	    tprime2[j] = 0.0;

	    uvprime[j] = 0.0;
	    vwprime[j] = 0.0;
	    uwprime[j] = 0.0;
	}

	FOR_Z{
	    FOR_Y{
		FOR_X{
		    int ii = GET3DINDEX_XYZ;
		    rhoAvg[j] += rho1[ii]/((double)(Nx*Nz));
		    uAvg[j]   += U[ii]/((double)(Nx*Nz));
		    vAvg[j]   += V[ii]/((double)(Nx*Nz));
		    wAvg[j]   += W[ii]/((double)(Nx*Nz));
		    pAvg[j]   += p[ii]/((double)(Nx*Nz));
		    tAvg[j]   += T[ii]/((double)(Nx*Nz));
		}
	    }
	}

	FOR_Z{
	    FOR_Y{
	        FOR_X{
		    int ii = GET3DINDEX_XYZ;
		    rhoprime2[j] += (rho1[ii]-rhoAvg[j])*(rho1[ii]-rhoAvg[j])/((double)(Nx*Nz)); 
		    uprime2[j] += (U[ii]-uAvg[j])*(U[ii]-uAvg[j])/((double)(Nx*Nz)); 
		    vprime2[j] += (V[ii]-vAvg[j])*(V[ii]-vAvg[j])/((double)(Nx*Nz)); 
		    wprime2[j] += (W[ii]-wAvg[j])*(W[ii]-wAvg[j])/((double)(Nx*Nz));
		    pprime2[j] += (p[ii]-pAvg[j])*(p[ii]-pAvg[j])/((double)(Nx*Nz));
		    tprime2[j] += (T[ii]-tAvg[j])*(T[ii]-tAvg[j])/((double)(Nx*Nz));
		    uvprime[j] += (U[ii]-uAvg[j])*(V[ii]-vAvg[j])/((double)(Nx*Nz)); 
		    uwprime[j] += (U[ii]-uAvg[j])*(W[ii]-wAvg[j])/((double)(Nx*Nz)); 
		    uwprime[j] += (V[ii]-vAvg[j])*(W[ii]-wAvg[j])/((double)(Nx*Nz)); 
	        }
	    }
	}
	
	
        ofstream outfile;
        outfile.precision(17);
        string outputFileName;

        outputFileName = "rhoprime2.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
	outfile << time << " ";
	FOR_Y{
            outfile << rhoprime2[j] << " ";
	}
	outfile << endl;
        outfile.close();

	
        outputFileName = "uprime2.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
	outfile << time << " ";
	FOR_Y{
            outfile << uprime2[j] << " ";
	}
	outfile << endl;
        outfile.close();


        outputFileName = "vprime2.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
	outfile << time << " ";
	FOR_Y{
            outfile << vprime2[j] << " ";
	}
	outfile << endl;
        outfile.close();


        outputFileName = "wprime2.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
	outfile << time << " ";
	FOR_Y{
            outfile << wprime2[j] << " ";
	}
	outfile << endl;
        outfile.close();


        outputFileName = "pprime2.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
	outfile << time << " ";
	FOR_Y{
            outfile << pprime2[j] << " ";
	}
	outfile << endl;
        outfile.close();


        outputFileName = "tprime2.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
	outfile << time << " ";
	FOR_Y{
            outfile << tprime2[j] << " ";
	}
	outfile << endl;
        outfile.close();


        outputFileName = "uvprime.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
	outfile << time << " ";
	FOR_Y{
            outfile << uvprime[j] << " ";
	}
	outfile << endl;
        outfile.close();

        outputFileName = "uwprime.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
	outfile << time << " ";
	FOR_Y{
            outfile << uwprime[j] << " ";
	}
	outfile << endl;
        outfile.close();

        outputFileName = "vwprime.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
	outfile << time << " ";
	FOR_Y{
            outfile << vwprime[j] << " ";
	}
	outfile << endl;
        outfile.close();

        outputFileName = "uAvg.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
	outfile << time << " ";
	FOR_Y{
            outfile << uAvg[j] << " ";
	}
	outfile << endl;
        outfile.close();

        outputFileName = "rhoAvg.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
	outfile << time << " ";
	FOR_Y{
            outfile << rhoAvg[j] << " ";
	}
	outfile << endl;
        outfile.close();



    }


    if(useTiming){
        cout << " > shrLyrInfo Timing: ";
        toc();
    }



};

*/
void UniformCSolver::checkSolution(){


    if(timeStep%ts->checkStep == 0){

	t2 = MPI_Wtime();

	IF_RANK0{
	    cout << endl;
            cout << "-------------------------------------------------" << endl;
            cout << " Step = "<< timeStep << ", time = " << time << ", dt = " << ts->dt << endl;
            cout << "-------------------------------------------------" << endl;
            cout << "  Time since last timestep = " << t2 - t1  << endl;
	}
	
        getRange(rho1, "RHO", pxSize[0], pxSize[1], pxSize[2], mpiRank);
        getRange(U, "U", pxSize[0], pxSize[1], pxSize[2], mpiRank);
        getRange(V, "V", pxSize[0], pxSize[1], pxSize[2], mpiRank);
        getRange(W, "W", pxSize[0], pxSize[1], pxSize[2], mpiRank);
        getRange(p, "P", pxSize[0], pxSize[1], pxSize[2], mpiRank);
        getRange(T, "T", pxSize[0], pxSize[1], pxSize[2], mpiRank);
        getRange(mu, "mu", pxSize[0], pxSize[1], pxSize[2], mpiRank);
        getRange(rhoE1, "RHOE", pxSize[0], pxSize[1], pxSize[2], mpiRank);
        getRange(sos, "SOS", pxSize[0], pxSize[1], pxSize[2], mpiRank);
        cout << endl;

	t1 = MPI_Wtime();
    }


};


void UniformCSolver::dumpSolution(){


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
	outputFileName.append(to_string(timeStep));

	MPI_File fh;
	MPI_Offset disp, filesize;

	MPI_File_open(MPI_COMM_WORLD, outputFileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	filesize = 0;
	MPI_File_set_size(fh, filesize);

	disp = 0;
	c2d->writeVar(fh, disp, 0, rho1);
	c2d->writeVar(fh, disp, 0, rhoU1);
	c2d->writeVar(fh, disp, 0, rhoV1);
	c2d->writeVar(fh, disp, 0, rhoW1);
	c2d->writeVar(fh, disp, 0, rhoE1);

	MPI_File_close(&fh);

	double time2 = MPI_Wtime();

	IF_RANK0{
	    cout << endl;
	    cout << " > File dump took " << time2-time1 << " seconds." << endl;
	}

}
/*
void UniformCSolver::writeImages(){

    if(useTiming) tic();

    if(timeStep%25==0){
	cout << " > Dumping images..." << endl;

	string timeStepString = to_string(timeStep);
	int zeroPad = 6;
	timeStepString = string(zeroPad - timeStepString.length(), '0') + timeStepString;

	//going to do our images in greyscale
	double dataMin, dataMax;
	
	//Density...
	getRangeValue(rho1, Nx, Ny, Nz, dataMin, dataMax);

	FOR_Y{
	    FOR_X{
		int k = (int)(Nz/2.0);
		int ii = GET3DINDEX_XYZ; 
		double f = (rho1[ii] - dataMin)/(dataMax - dataMin);
		int g = (int)(f*255.0);
		pngXY->set(i,j,g,g,g);
	    }
	}
	string imageName = "imagesXY/rhoXY.";
	imageName.append(timeStepString);
	imageName.append(".png");
	pngXY->write(imageName.c_str()); 

	//Pressure...
	getRangeValue(p, Nx, Ny, Nz, dataMin, dataMax);

	FOR_Y{
	    FOR_X{
		int k = (int)(Nz/2.0);
		int ii = GET3DINDEX_XYZ; 
		double f = (p[ii] - dataMin)/(dataMax - dataMin);
		int g = (int)(f*255.0);
		pngXY->set(i,j,g,g,g);
	    }
	}
	imageName = "imagesXY/pXY.";
	imageName.append(timeStepString);
	imageName.append(".png");
	pngXY->write(imageName.c_str()); 

	//U...
	getRangeValue(U, Nx, Ny, Nz, dataMin, dataMax);

	FOR_Y{
	    FOR_X{
		int k = (int)(Nz/2.0);
		int ii = GET3DINDEX_XYZ; 
		double f = (U[ii] - dataMin)/(dataMax - dataMin);
		int g = (int)(f*255.0);
		pngXY->set(i,j,g,g,g);
	    }
	}
	imageName = "imagesXY/UXY.";
	imageName.append(timeStepString);
	imageName.append(".png");
	pngXY->write(imageName.c_str()); 

	//V...
	getRangeValue(V, Nx, Ny, Nz, dataMin, dataMax);

	FOR_Y{
	    FOR_X{
		int k = (int)(Nz/2.0);
		int ii = GET3DINDEX_XYZ; 
		double f = (V[ii] - dataMin)/(dataMax - dataMin);
		int g = (int)(f*255.0);
		pngXY->set(i,j,g,g,g);
	    }
	}
	imageName = "imagesXY/VXY.";
	imageName.append(timeStepString);
	imageName.append(".png");
	pngXY->write(imageName.c_str()); 

	//Dilataion...
	double *dil = new double[N];

	#pragma omp parallel for
	FOR_XYZ dil[ip] = Ux[ip] + Vy[ip] + Wz[ip];

	getRangeValue(dil, Nx, Ny, Nz, dataMin, dataMax);

	FOR_Y{
	    FOR_X{
		int k = (int)(Nz/2.0);
		int ii = GET3DINDEX_XYZ; 
		double f = (dil[ii] - dataMin)/(dataMax - dataMin);
		int g = (int)(f*255.0);
		pngXY->set(i,j,g,g,g);
	    }
	}
	imageName = "imagesXY/dilXY.";
	imageName.append(timeStepString);
	imageName.append(".png");
	pngXY->write(imageName.c_str()); 

	delete[] dil;

	//XZ Plane,
	//U
	getRangeValue(U, Nx, Ny, Nz, dataMin, dataMax);

	FOR_X{
	    FOR_Z{
		int j = (int)(Ny/2.0);
		int ii = GET3DINDEX_XYZ; 
		double f = (U[ii] - dataMin)/(dataMax - dataMin);
		int g = (int)(f*255.0);
		pngXZ->set(k,i,g,g,g);
	    }
	}
	imageName = "imagesXZ/UXZ.";
	imageName.append(timeStepString);
	imageName.append(".png");
	pngXZ->write(imageName.c_str()); 

	//V
	getRangeValue(V, Nx, Ny, Nz, dataMin, dataMax);

	FOR_X{
	    FOR_Z{
		int j = (int)(Ny/2.0);
		int ii = GET3DINDEX_XYZ; 
		double f = (V[ii] - dataMin)/(dataMax - dataMin);
		int g = (int)(f*255.0);
		pngXZ->set(k,i,g,g,g);
	    }
	}
	imageName = "imagesXZ/VXZ.";
	imageName.append(timeStepString);
	imageName.append(".png");
	pngXZ->write(imageName.c_str()); 

	//P
	getRangeValue(p, Nx, Ny, Nz, dataMin, dataMax);

	FOR_X{
	    FOR_Z{
		int j = (int)(Ny/2.0);
		int ii = GET3DINDEX_XYZ; 
		double f = (p[ii] - dataMin)/(dataMax - dataMin);
		int g = (int)(f*255.0);
		pngXZ->set(k,i,g,g,g);
	    }
	}
	imageName = "imagesXZ/pXZ.";
	imageName.append(timeStepString);
	imageName.append(".png");
	pngXZ->write(imageName.c_str()); 



    }

    if(useTiming){
        cout << " > writeImage Timing: ";
        toc();
    }


}
*/
void UniformCSolver::checkEnd(){


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
/*
void UniformCSolver::reportAll(){

   cout << "REPORT ALL" << endl;

   getRange(Ux, "Ux", Nx, Ny, Nz);
   getRange(Uxx, "Uxx", Nx, Ny, Nz);
   getRange(Uy, "Uy", Nx, Ny, Nz);
   getRange(Uyy, "Uyy", Nx, Ny, Nz);
   getRange(Uz, "Uz", Nx, Ny, Nz);
   getRange(Uzz, "Uzz", Nx, Ny, Nz);
   getRange(Uxy, "Uxy", Nx, Ny, Nz);
   getRange(Uyz, "Uyz", Nx, Ny, Nz);
   getRange(Uxz, "Uxz", Nx, Ny, Nz);
cout << " " << endl;
   getRange(Vx, "Vx", Nx, Ny, Nz);
   getRange(Vxx, "Vxx", Nx, Ny, Nz);
   getRange(Vy, "Vy", Nx, Ny, Nz);
   getRange(Vyy, "Vyy", Nx, Ny, Nz);
   getRange(Vz, "Vz", Nx, Ny, Nz);
   getRange(Vzz, "Vzz", Nx, Ny, Nz);
   getRange(Vxy, "Vxy", Nx, Ny, Nz);
   getRange(Vyz, "Vyz", Nx, Ny, Nz);
   getRange(Vxz, "Vxz", Nx, Ny, Nz);
cout << " " << endl;
   getRange(Wx, "Wx", Nx, Ny, Nz);
   getRange(Wxx, "Wxx", Nx, Ny, Nz);
   getRange(Wy, "Wy", Nx, Ny, Nz);
   getRange(Wyy, "Wyy", Nx, Ny, Nz);
   getRange(Wz, "Wz", Nx, Ny, Nz);
   getRange(Wzz, "Wzz", Nx, Ny, Nz);
   getRange(Wxy, "Wxy", Nx, Ny, Nz);
   getRange(Wyz, "Wyz", Nx, Ny, Nz);
   getRange(Wxz, "Wxz", Nx, Ny, Nz);
cout << " " << endl;
   getRange(Tx, "Ux", Nx, Ny, Nz);
   getRange(Txx, "Uxx", Nx, Ny, Nz);
   getRange(Ty, "Uy", Nx, Ny, Nz);
   getRange(Tyy, "Uyy", Nx, Ny, Nz);
   getRange(Tz, "Uz", Nx, Ny, Nz);
   getRange(Tzz, "Uzz", Nx, Ny, Nz);
cout << " " << endl;
   getRange(contEulerX, "contEulerX", Nx, Ny, Nz);
   getRange(contEulerY, "contEulerY", Nx, Ny, Nz);
   getRange(contEulerZ, "contEulerZ", Nx, Ny, Nz);
cout << " " << endl;
   getRange(momXEulerX, "momXEulerX", Nx, Ny, Nz);
   getRange(momXEulerY, "momXEulerY", Nx, Ny, Nz);
   getRange(momXEulerZ, "momXEulerZ", Nx, Ny, Nz);
cout << " " << endl;
   getRange(momYEulerX, "momYEulerX", Nx, Ny, Nz);
   getRange(momYEulerY, "momYEulerY", Nx, Ny, Nz);
   getRange(momYEulerZ, "momYEulerZ", Nx, Ny, Nz);
cout << " " << endl;
   getRange(momZEulerX, "momZEulerX", Nx, Ny, Nz);
   getRange(momZEulerY, "momZEulerY", Nx, Ny, Nz);
   getRange(momZEulerZ, "momZEulerZ", Nx, Ny, Nz);
cout << " " << endl;
   getRange(engyEulerX, "engyEulerX", Nx, Ny, Nz);
   getRange(engyEulerY, "engyEulerY", Nx, Ny, Nz);
   getRange(engyEulerZ, "engyEulerZ", Nx, Ny, Nz);
cout << " " << endl;
   getRange(rho1, "rho1", Nx, Ny, Nz);
   getRange(rhok, "rhok", Nx, Ny, Nz);
   getRange(rhok2, "rhok2", Nx, Ny, Nz);
   getRange(rho2, "rho2", Nx, Ny, Nz);
cout << " " << endl;
   getRange(rhoU1, "rhoU1", Nx, Ny, Nz);
   getRange(rhoUk, "rhoUk", Nx, Ny, Nz);
   getRange(rhoUk2, "rhoUk2", Nx, Ny, Nz);
   getRange(rhoU2, "rhoU2", Nx, Ny, Nz);
cout << " " << endl;
   getRange(rhoV1, "rhoV1", Nx, Ny, Nz);
   getRange(rhoVk, "rhoVk", Nx, Ny, Nz);
   getRange(rhoVk2, "rhoVk2", Nx, Ny, Nz);
   getRange(rhoV2, "rhoV2", Nx, Ny, Nz);
cout << " " << endl;
   getRange(rhoW1, "rhoW1", Nx, Ny, Nz);
   getRange(rhoWk, "rhoWk", Nx, Ny, Nz);
   getRange(rhoWk2, "rhoWk2", Nx, Ny, Nz);
   getRange(rhoW2, "rhoW2", Nx, Ny, Nz);
cout << " " << endl;
   getRange(rhoE1, "rhoE1", Nx, Ny, Nz);
   getRange(rhoEk, "rhoEk", Nx, Ny, Nz);
   getRange(rhoEk2, "rhoEk2", Nx, Ny, Nz);
   getRange(rhoE2, "rhoE2", Nx, Ny, Nz);
cout << " " << endl;
   getRange(p, "p", Nx, Ny, Nz);
   getRange(U, "U", Nx, Ny, Nz);
   getRange(V, "V", Nx, Ny, Nz);
   getRange(W, "W", Nx, Ny, Nz);
   getRange(T, "T", Nx, Ny, Nz);
   getRange(mu, "mu", Nx, Ny, Nz);
   getRange(Amu, "Amu", Nx, Ny, Nz);
   getRange(sos, "sos", Nx, Ny, Nz);
cout << " " << endl;

}
*/
/////////////////////////////////////
//Our Generalized Solver Functions //
/////////////////////////////////////

void UniformCSolver::preStep(){
 
   if(timeStep == 0){
        dumpSolution();
//	writeImages();
    }
    
    calcDtFromCFL();

}

void UniformCSolver::preSubStep(){

    preStepBCHandling();
    preStepDerivatives();

}

void UniformCSolver::solveEqnSet(){

    solveContinuity();
    solveXMomentum();
    solveYMomentum();
    solveZMomentum();
    solveEnergy();

}

void UniformCSolver::postSubStep(){

    postStepBCHandling();

}

void UniformCSolver::updateData(){

    if(rkLast){
        filterConservedData();
    }
    updateNonConservedData();

}

void UniformCSolver::postStep(){

    //calcTurbulenceQuantities();
    //shearLayerInfoCalc();
    //calcTaylorGreenQuantities();

    updateSponge();
    checkSolution();

    if(timeStep%ts->dumpStep == 0)
        dumpSolution();

    checkEnd();
    //reportAll();

}

