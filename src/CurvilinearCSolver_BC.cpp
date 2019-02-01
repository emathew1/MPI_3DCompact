#include "CurvilinearCSolver.hpp"

void CurvilinearCSolver::setInitialConditions(){


    if(useTiming) ft1 = MPI_Wtime();

    IF_RANK0{
        cout << endl;
        cout << " > Setting initial conditions..." << endl; 
    }


    //Initialize sponge boundary conditions if necessary
    spg = NULL;

    //TODO NEED TO ADD EXCEPTIONS FOR WHEN WRONG KIND OF SPONGE OR UNIMPLEMENTED SPONGE BC IS TRYING TO BE USED
    if( bc->bcX0 == Options::SPONGE || \
        bc->bcX1 == Options::SPONGE || \
        bc->bcY0 == Options::SPONGE || \
        bc->bcY1 == Options::SPONGE || \
        bc->bcZ0 == Options::SPONGE || \
        bc->bcZ1 == Options::SPONGE ){
        spongeFlag = true;
        spg = new SpongeBC(msh, dom, ig, bc, c2d, opt, mpiRank);
    }else{
        spongeFlag = false;
    }

    //just do the simple stuff in a loop...
    FOR_XYZ_YPEN{
	U[ip]	 = U0[ip];
	V[ip] 	 = V0[ip];
	W[ip] 	 = W0[ip];
	p[ip]	 = p0[ip];

	rho1[ip]  = rho0[ip];
	rhoU1[ip] = rho1[ip]*U[ip];	
	rhoV1[ip] = rho1[ip]*V[ip];	
	rhoW1[ip] = rho1[ip]*W[ip];

	Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 
    }

    //Can now release initial condition data...
    c2d->deallocXYZ(rho0);
    c2d->deallocXYZ(U0);
    c2d->deallocXYZ(V0);
    c2d->deallocXYZ(W0);
    c2d->deallocXYZ(p0);

    //Call the ideal gas relations for the slightly more involved stuff..
    FOR_XYZ_YPEN{
	rhoE1[ip] = ig->solverhoE(rho1[ip], p[ip], U[ip], V[ip], W[ip]);
    	T[ip]     = ig->solveT(rho1[ip], p[ip]);
        mu[ip]    = ig->solveMu(T[ip]);
        sos[ip]   = ig->solveSOS(rho1[ip], p[ip]);
    }
    //This is where we'll do the boundary condition specific stuff...
    bool wallBCFlag = false;

    //--------------------------------
    //No-Slip Wall Boundary Conditions
    //--------------------------------

    if(bc->bcX0 == Options::ADIABATIC_WALL || bc->bcX0 == Options::CONST_T_WALL){
	wallBCFlag = true;
        FOR_X0_YPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
	    if(bc->bcX0 == Options::ADIABATIC_WALL){

		int index[10] = FILL_GETMAJIND_YPEN_Xp;
		double T_out[10];
		getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivX->calcNeumann(T_out);

	    }
        }END_FORX0
    }


    if(bc->bcX1 == Options::ADIABATIC_WALL || bc->bcX1 == Options::CONST_T_WALL){
	wallBCFlag = true;
        FOR_X1_YPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
	    if(bc->bcX1 == Options::CONST_T_WALL){

		int index[10] = FILL_GETMAJIND_YPEN_Xm;
	 	double T_out[10];
		getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivX->calcNeumann(T_out);

	    }
        }END_FORX1
    }

    if(bc->bcY0 == Options::ADIABATIC_WALL || bc->bcY0 == Options::CONST_T_WALL){
	wallBCFlag = true;
        FOR_Y0_YPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
	    if(bc->bcY0 == Options::ADIABATIC_WALL){

		int index[10] = FILL_GETMAJIND_YPEN_Yp;
		double T_out[10];
		getDataFromIndex(T, index, 10, T_out);
	        T[ip] = derivY->calcNeumann(T_out);

	    }
        }END_FORY0
    }

    if(bc->bcY1 == Options::ADIABATIC_WALL || bc->bcY1 == Options::CONST_T_WALL){
	wallBCFlag = true;
        FOR_Y1_YPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
	    if(bc->bcY1 == Options::ADIABATIC_WALL){

		int index[10] = FILL_GETMAJIND_YPEN_Ym;
		double T_out[10];
		getDataFromIndex(T, index, 10, T_out);
	        T[ip] = derivY->calcNeumann(T_out);

	    }
        }END_FORY1
    }

    if(bc->bcZ0 == Options::ADIABATIC_WALL || bc->bcZ0 == Options::CONST_T_WALL){
	wallBCFlag = true;
        FOR_Z0_YPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
	    if(bc->bcZ0 == Options::ADIABATIC_WALL){

		int index[10] = FILL_GETMAJIND_YPEN_Zp;
		double T_out[10];
		getDataFromIndex(T, index, 10, T_out);
		T[ip] = derivZ->calcNeumann(T_out);

	    }
        }END_FORZ0

    }

    if(bc->bcZ1 == Options::ADIABATIC_WALL || bc->bcZ1 == Options::CONST_T_WALL){
	wallBCFlag = true;
        FOR_Z1_YPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
	    if(bc->bcZ1 == Options::ADIABATIC_WALL){
		
		int index[10] = FILL_GETMAJIND_YPEN_Zm;
		double T_out[10];
		getDataFromIndex(T, index, 10, T_out);
	        T[ip] = derivZ->calcNeumann(T_out);
	    }

        }END_FORZ1

    }

    //-------------------------------
    //Moving Wall Boundary Conditions
    //-------------------------------

    if(bc->bcX0 == Options::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_X0_YPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = X0WallV;
            W[ip]  = X0WallW;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = rho1[ip]*X0WallV;
            rhoW1[ip] = rho1[ip]*X0WallW;

	    int index[10] = FILL_GETMAJIND_YPEN_Xp;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivX->calcNeumann(T_out);

        }END_FORX0
    }

    if(bc->bcX1 == Options::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_X1_YPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = X1WallV;
            W[ip]  = X1WallW;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = rho1[ip]*X1WallV;
            rhoW1[ip] = rho1[ip]*X1WallW;

	    int index[10] = FILL_GETMAJIND_YPEN_Xm;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivX->calcNeumann(T_out);

        }END_FORX1
    }

    if(bc->bcY0 == Options::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Y0_YPEN_MAJ{
            U[ip]  = Y0WallU;
            V[ip]  = 0.0;
            W[ip]  = Y0WallW;
            rhoU1[ip] = rho1[ip]*Y0WallU;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = rho1[ip]*Y0WallW;

	    int index[10] = FILL_GETMAJIND_YPEN_Yp;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivY->calcNeumann(T_out);

        }END_FORY0

   }

    if(bc->bcY1 == Options::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Y1_YPEN_MAJ{
            U[ip]  = Y1WallU;
            V[ip]  = 0.0;
            W[ip]  = Y1WallW;
            rhoU1[ip] = rho1[ip]*Y1WallU;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = rho1[ip]*Y1WallW;

	    int index[10] = FILL_GETMAJIND_YPEN_Ym;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivY->calcNeumann(T_out);

        }END_FORY1

   }

    if(bc->bcZ0 == Options::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Z0_YPEN_MAJ{
            U[ip]  = Z0WallU;
            V[ip]  = Z0WallV;
            W[ip]  = 0.0;
            rhoU1[ip] = rho1[ip]*Z0WallU;
            rhoV1[ip] = rho1[ip]*Z0WallV;
            rhoW1[ip] = 0.0;

	    int index[10] = FILL_GETMAJIND_YPEN_Zp;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivZ->calcNeumann(T_out);

        }END_FORZ0


    }

    if(bc->bcZ1 == Options::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Z1_YPEN_MAJ{
            U[ip]  = Z1WallU;
            V[ip]  = Z1WallV;
            W[ip]  = 0.0;
            rhoU1[ip] = rho1[ip]*Z1WallU;
            rhoV1[ip] = rho1[ip]*Z1WallV;
            rhoW1[ip] = 0.0;

	    int index[10] = FILL_GETMAJIND_YPEN_Zm;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivZ->calcNeumann(T_out);

        }END_FORZ1

    }




    if(wallBCFlag == true){
	//Need to update the pressure, sos, and rhoE fields at the boundaries with walls...
	FOR_XYZ_YPEN{
	    p[ip]     = ig->solvep_idealgas(rho1[ip], T[ip]);
	    sos[ip]   = ig->solveSOS(rho1[ip],p[ip]);
	    rhoE1[ip] = ig->solverhoE(rho1[ip], p[ip], U[ip], V[ip], W[ip]);

	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	    Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 
	}
    }

    if(spongeFlag == true){

        if(opt->spongeFromRestart){

            spg->readSpongeAvgFromRestart();   

    	    getRange(spg->spongeRhoAvg,   "spongeRhoAvg",  pySize[0], pySize[1], pySize[2], mpiRank);
    	    getRange(spg->spongeRhoUAvg,  "spongeRhoUAvg", pySize[0], pySize[1], pySize[2], mpiRank);
    	    getRange(spg->spongeRhoVAvg,  "spongeRhoVAvg", pySize[0], pySize[1], pySize[2], mpiRank);
    	    getRange(spg->spongeRhoWAvg,  "spongeRhoWAvg", pySize[0], pySize[1], pySize[2], mpiRank);
    	    getRange(spg->spongeRhoEAvg,  "spongeRhoEAvg", pySize[0], pySize[1], pySize[2], mpiRank);

	}else{
	    FOR_XYZ_YPEN{
	        spg->spongeRhoAvg[ip]  = rho1[ip];
	        spg->spongeRhoUAvg[ip] = rhoU1[ip];
	        spg->spongeRhoVAvg[ip] = rhoV1[ip];
	        spg->spongeRhoWAvg[ip] = rhoW1[ip];
	        spg->spongeRhoEAvg[ip] = rhoE1[ip];
 	    }
	}
    }



    IF_RANK0 std::cout << " > Finished initialization of flow field " << std::endl;

    getRange(rho1,  "rho0", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(rhoU1, "rhoU0", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(rhoV1, "rhoV0", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(rhoW1, "rhoW0", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(rhoE1, "rhoE0", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(T, "T0", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(p, "p0", pySize[0], pySize[1], pySize[2], mpiRank);

    //Run the initial hook
    initialHook();

    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > setInitCond Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }


}



void CurvilinearCSolver::preStepBCHandling(){

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

    if(bc->bcX0 == Options::ADIABATIC_WALL || bc->bcX0 == Options::CONST_T_WALL){
	FOR_X0_YPEN_MAJ{
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    if(bc->bcX0 == Options::ADIABATIC_WALL){
	        int index[10] = FILL_GETMAJIND_YPEN_Xp;
	        double T_out[10];
	        getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivX->calcNeumann(T_out);
	    }
	   p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	   rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
  	   Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	   Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	   Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 

	}END_FORX0

    }
    

    if(bc->bcX1 == Options::ADIABATIC_WALL || bc->bcX1 == Options::CONST_T_WALL){
	FOR_X1_YPEN_MAJ{
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    if(bc->bcX1 == Options::ADIABATIC_WALL){
	        int index[10] = FILL_GETMAJIND_YPEN_Xm;
	        double T_out[10];
	        getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivX->calcNeumann(T_out);
	    }
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	    Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 

	}END_FORX1
    }   


    if(bc->bcY0 == Options::ADIABATIC_WALL || bc->bcY0 == Options::CONST_T_WALL){

	FOR_Y0_YPEN_MAJ{
	    if(bc->bcY0 == Options::ADIABATIC_WALL){
	        int index[10] = FILL_GETMAJIND_YPEN_Yp;
	        double T_out[10];
	        getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivY->calcNeumann(T_out);
	    }
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]);
   	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	    Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 
	}END_FORY0
    }
    

    if(bc->bcY1 == Options::ADIABATIC_WALL || bc->bcY1 == Options::CONST_T_WALL){

	FOR_Y1_YPEN_MAJ{

	    if(bc->bcY1 == Options::ADIABATIC_WALL){
	        int index[10] = FILL_GETMAJIND_YPEN_Ym;
	        double T_out[10];
	        getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivY->calcNeumann(T_out);
	    }
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	    Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 

	}END_FORY1
    }



    if(bc->bcZ0 == Options::ADIABATIC_WALL || bc->bcZ0 == Options::CONST_T_WALL){

	FOR_Z0_YPEN_MAJ{
	    if(bc->bcZ0 == Options::ADIABATIC_WALL){
	        int index[10] = FILL_GETMAJIND_YPEN_Zp;
	        double T_out[10];
	        getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivZ->calcNeumann(T_out);
	    }
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	    Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 
	}END_FORZ0
    }



    if(bc->bcZ1 == Options::ADIABATIC_WALL || bc->bcZ1 == Options::CONST_T_WALL){
	FOR_Z1_YPEN_MAJ{
	    if(bc->bcZ1 == Options::ADIABATIC_WALL){
	        int index[10] = FILL_GETMAJIND_YPEN_Zm;
	        double T_out[10];
	        getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivZ->calcNeumann(T_out);
	    }
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	    Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 
	}END_FORZ1
    }


    //-------------------------------
    //Moving wall boundary conditions
    //-------------------------------

    if(bc->bcX0 == Options::MOVING_ADIABATIC_WALL){
        FOR_X0_YPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = X0WallV;
            W[ip]  = X0WallW;
            rhoUP[ip] = 0.0;
            rhoVP[ip] = rhoP[ip]*X0WallV;
            rhoWP[ip] = rhoP[ip]*X0WallW;

	    int index[10] = FILL_GETMAJIND_YPEN_Xp;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivX->calcNeumann(T_out);

	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	    Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 
        }END_FORX0
    }

    if(bc->bcX1 == Options::MOVING_ADIABATIC_WALL){
        FOR_X1_YPEN_MAJ{
            U[ip]  = 0.0;
            V[ip]  = X1WallV;
            W[ip]  = X1WallW;
            rhoUP[ip] = 0.0;
            rhoVP[ip] = rhoP[ip]*X1WallV;
            rhoWP[ip] = rhoP[ip]*X1WallW;

	    int index[10] = FILL_GETMAJIND_YPEN_Xm;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivX->calcNeumann(T_out);

	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	    Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 
        }END_FORX1
    }

    if(bc->bcY0 == Options::MOVING_ADIABATIC_WALL){

        FOR_Y0_YPEN_MAJ{

	    int index[10] = FILL_GETMAJIND_YPEN_Yp;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivY->calcNeumann(T_out);

            U[ip]  = Y0WallU;
            V[ip]  = 0.0;
            W[ip]  = Y0WallW;
            rhoUP[ip] = rhoP[ip]*Y0WallU;
            rhoVP[ip] = 0.0;
            rhoWP[ip] = rhoP[ip]*Y0WallW;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	    Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 
        }END_FORY0
    }

    if(bc->bcY1 == Options::MOVING_ADIABATIC_WALL){

        FOR_Y1_YPEN_MAJ{

	    int index[10] = FILL_GETMAJIND_YPEN_Ym;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivY->calcNeumann(T_out);

            U[ip]  = Y1WallU;
            V[ip]  = 0.0;
            W[ip]  = Y1WallW;
            rhoUP[ip] = rhoP[ip]*Y1WallU;
            rhoVP[ip] = 0.0;
            rhoWP[ip] = rhoP[ip]*Y1WallW;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	    Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 
        }END_FORY1
    }

    if(bc->bcZ0 == Options::MOVING_ADIABATIC_WALL){

        FOR_Z0_YPEN_MAJ{
	
	    int index[10] = FILL_GETMAJIND_YPEN_Zp;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivZ->calcNeumann(T_out);

            U[ip]  = Z0WallU;
            V[ip]  = Z0WallV;
            W[ip]  = 0.0;
            rhoUP[ip] = rhoP[ip]*Z0WallU;
            rhoVP[ip] = rhoP[ip]*Z0WallV;
            rhoWP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	    Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 
        }END_FORZ0
    }

    if(bc->bcZ1 == Options::MOVING_ADIABATIC_WALL){

        FOR_Z1_YPEN_MAJ{
	    
	    int index[10] = FILL_GETMAJIND_YPEN_Zm;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivZ->calcNeumann(T_out);

            U[ip]  = Z1WallU;
            V[ip]  = Z1WallV;
            W[ip]  = 0.0;
            rhoUP[ip] = rhoP[ip]*Z1WallU;
            rhoVP[ip] = rhoP[ip]*Z1WallV;
            rhoWP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip] + J13[ip]*W[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip] + J23[ip]*W[ip]; 
	    Wcurv[ip] = J31[ip]*U[ip] + J32[ip]*V[ip] + J33[ip]*W[ip]; 
        }END_FORZ1
    }


    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > preBCHandl Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }



}



void CurvilinearCSolver::postStepBCHandling(){

    if(useTiming) ft1 = MPI_Wtime();

    ////////////////////////////////
    //ADIABATIC AND MOVING WALL BC// 
    ////////////////////////////////

    if(bc->bcX0 == Options::ADIABATIC_WALL || bc->bcX0 == Options::MOVING_ADIABATIC_WALL || bc->bcX0 == Options::CONST_T_WALL){
	FOR_X0_YPEN_MAJ{
	    rhok2[ip]  = -ts->dt*J[ip]*cont_1[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORX0
    }

    if(bc->bcX1 == Options::ADIABATIC_WALL || bc->bcX1 == Options::MOVING_ADIABATIC_WALL || bc->bcX1 == Options::CONST_T_WALL){
	FOR_X1_YPEN_MAJ{
	    rhok2[ip]  = -ts->dt*J[ip]*cont_1[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORX1
    }   

    if(bc->bcY0 == Options::ADIABATIC_WALL || bc->bcY0 == Options::MOVING_ADIABATIC_WALL || bc->bcY0 == Options::CONST_T_WALL){
	FOR_Y0_YPEN_MAJ{
	    rhok2[ip]  = -ts->dt*J[ip]*cont_2[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORY0
    }

    if(bc->bcY1 == Options::ADIABATIC_WALL || bc->bcY1 == Options::MOVING_ADIABATIC_WALL || bc->bcY1 == Options::CONST_T_WALL){
	FOR_Y1_YPEN_MAJ{
	    rhok2[ip]  = -ts->dt*J[ip]*cont_2[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORY1
    }

    if(bc->bcZ0 == Options::ADIABATIC_WALL || bc->bcZ0 == Options::MOVING_ADIABATIC_WALL || bc->bcZ0 == Options::CONST_T_WALL){
	FOR_Z0_YPEN_MAJ{
	    rhok2[ip]  = -ts->dt*J[ip]*cont_3[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORZ0
    }

    if(bc->bcZ1 == Options::ADIABATIC_WALL || bc->bcZ1 == Options::MOVING_ADIABATIC_WALL || bc->bcZ1 == Options::CONST_T_WALL){
	FOR_Z1_YPEN_MAJ{
	    rhok2[ip]  = -ts->dt*J[ip]*cont_3[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORZ1
    }


    /////////////
    //SPONGE BC//
    /////////////

    if(bc->bcX0 == Options::SPONGE){
	FOR_X0_YPEN_MAJ{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORX0
    }

    if(bc->bcX1 == Options::SPONGE){
	FOR_X1_YPEN_MAJ{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORX1
    }   

    if(bc->bcY0 == Options::SPONGE){
	FOR_Y0_YPEN_MAJ{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORY0
    }

    if(bc->bcY1 == Options::SPONGE){
	FOR_Y1_YPEN_MAJ{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORY1
    }

    if(bc->bcZ0 == Options::SPONGE){
	FOR_Z0_YPEN_MAJ{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORZ0
    }

    if(bc->bcZ1 == Options::SPONGE){
	FOR_Z1_YPEN_MAJ{
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


void CurvilinearCSolver::updateSponge(){

    if(useTiming) ft1 = MPI_Wtime();

    if(spongeFlag){
	double eps = 1.0/(spg->avgT/ts->dt + 1.0);
	FOR_XYZ_YPEN{
	    spg->spongeRhoAvg[ip]  += eps*(rho1[ip]  - spg->spongeRhoAvg[ip]);	
	    spg->spongeRhoUAvg[ip] += eps*(rhoU1[ip] - spg->spongeRhoUAvg[ip]);	
	    spg->spongeRhoVAvg[ip] += eps*(rhoV1[ip] - spg->spongeRhoVAvg[ip]);	
	    spg->spongeRhoWAvg[ip] += eps*(rhoW1[ip] - spg->spongeRhoWAvg[ip]);	
	    spg->spongeRhoEAvg[ip] += eps*(rhoE1[ip] - spg->spongeRhoEAvg[ip]);	
	    spg->spongeRhoEAvg[ip] = spg->epsP*spg->spongeRhoEAvg[ip] + (1.0 -  spg->epsP)*(spg->spongeP/(ig->gamma-1.0) \
					 + 0.5*(spg->spongeRhoUAvg[ip]*spg->spongeRhoUAvg[ip] + spg->spongeRhoVAvg[ip]*spg->spongeRhoVAvg[ip] \
					 + spg->spongeRhoWAvg[ip]*spg->spongeRhoWAvg[ip])/spg->spongeRhoAvg[ip]);
	}
	
        if(bc->bcX0 == Options::SPONGE){
	    FOR_X0_YPEN_MAJ{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORX0
        }

        if(bc->bcX1 == Options::SPONGE){
	    FOR_X1_YPEN_MAJ{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORX1
        }   

        if(bc->bcY0 == Options::SPONGE){
	    FOR_Y0_YPEN_MAJ{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORY0
        }

        if(bc->bcY1 == Options::SPONGE){
	    FOR_Y1_YPEN_MAJ{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORY1
        }

        if(bc->bcZ0 == Options::SPONGE){
	    FOR_Z0_YPEN_MAJ{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORZ0

        }

        if(bc->bcZ1 == Options::SPONGE){
	    FOR_Z1_YPEN_MAJ{
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


