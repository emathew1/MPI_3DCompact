#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <fstream>
#include <assert.h>

using namespace std;

#include "C2Decomp.hpp"
#include "Macros.hpp"
#include "TimeStepping.hpp"
#include "Domain.hpp"
#include "BC.hpp" 

#include "AbstractCSolver.hpp"
#include "CurvilinearCSolver.hpp"

#include "AbstractSingleBlockMesh.hpp"
#include "AlgebraicSingleBlockMesh.hpp"

#include "AbstractRK.hpp"
#include "TVDRK3.hpp"

#include "CurvilinearInterpolator.hpp"

int main(int argc, char *argv[]){
   int ierr, totRank, mpiRank;

    //Initialize MPI
    ierr = MPI_Init( &argc, &argv);

    //Get the number of processes
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &totRank);

    //Get the local rank
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    IF_RANK0{
        cout << endl;
        cout << " --------------------" << endl;
    	cout << "  Starting up Solver " << endl;
    	cout << " --------------------" << endl;
	cout << endl;
	cout << " > Running on " << totRank << " cores!" << endl;
	cout << endl;
    }
 

    ////////////////////////////////////
    //Time Stepping info intialization//
    ////////////////////////////////////
    TimeStepping::TimeSteppingType timeSteppingType = TimeStepping::CONST_CFL;
    double CFL       = 0.6;
    int maxTimeStep  = 25000;
    double maxTime   = 3000.0;
    int filterStep   = 1;
    int checkStep    = 1;
    int dumpStep     = 2500;
    int imageStep    = 25;
    TimeStepping *ts = new TimeStepping(timeSteppingType, 
					CFL, 
					maxTimeStep, 
					maxTime, 
					filterStep, 
					checkStep, 
					dumpStep,
					imageStep, 
					mpiRank);


    ///////////////////////////
    //Boundary Condition Info//
    ///////////////////////////
/*
    BC::BCType bcXType = BC::PERIODIC_SOLVE;
    BC::BCType bcYType = BC::PERIODIC_SOLVE;
    BC::BCType bcZType = BC::PERIODIC_SOLVE;

    BC::BCKind bcX0 = BC::PERIODIC;
    BC::BCKind bcX1 = BC::PERIODIC;
    BC::BCKind bcY0 = BC::PERIODIC;
    BC::BCKind bcY1 = BC::PERIODIC;
    BC::BCKind bcZ0 = BC::PERIODIC;
    BC::BCKind bcZ1 = BC::PERIODIC;
*/

    BC::BCType bcXType = BC::DIRICHLET_SOLVE;
    BC::BCType bcYType = BC::DIRICHLET_SOLVE;
    BC::BCType bcZType = BC::PERIODIC_SOLVE;

    BC::BCKind bcX0 = BC::RECT_CURVILINEARSPONGE;
    BC::BCKind bcX1 = BC::RECT_CURVILINEARSPONGE;
    BC::BCKind bcY0 = BC::ADIABATIC_WALL;
    BC::BCKind bcY1 = BC::ADIABATIC_WALL;
    BC::BCKind bcZ0 = BC::PERIODIC;
    BC::BCKind bcZ1 = BC::PERIODIC;


    bool periodicBC[3];
    BC *bc = new BC(bcXType, bcX0, bcX1,
                    bcYType, bcY0, bcY1,
                    bcZType, bcZ0, bcZ1,
		    periodicBC, mpiRank);

    /////////////////////////
    //Initialize the Domain//
    /////////////////////////
    int    Nx = 256,
           Ny = 100,
           Nz = 32;

    //For curvilinear coordinates these should all correspond to the max xi, eta, and zeta values
    double Lx = 1.0,
           Ly = 1.0,
           Lz = 1.0;
    Domain *d = new Domain(bc, Nx, Ny, Nz, Lx, Ly, Lz, mpiRank);


    /////////////////////////////
    //Initializing Pencil Decomp//
    //////////////////////////////
 
    int pRow = 0, 
	pCol = 0;
    IF_RANK0 cout << endl << " > Initializing the pencil decomposition... " << endl;

    C2Decomp *c2d;
    c2d = new C2Decomp(Nx, Ny, Nz, pRow, pCol, periodicBC);
  
    IF_RANK0 cout << " > Successfully initialized! " << endl;
    IF_RANK0 cout << " > Handing some decomp info back to the domain object... " << endl;
    d->setPencilDecompInfo(c2d);

    //This is info needed to use the macros...
    int pxSize[3], pySize[3], pzSize[3];
    int pxStart[3], pyStart[3], pzStart[3];
    int pxEnd[3], pyEnd[3], pzEnd[3];
    d->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);

    /////////////////////////
    //Initialize the Solver//
    /////////////////////////
    double alphaF  = 0.35;
    double mu_ref  = 0.00005;
    bool useTiming = false;
    AbstractCSolver *cs;
    cs = new CurvilinearCSolver(c2d, d, bc, ts, alphaF, mu_ref, useTiming);
    
    cs->msh = new AlgebraicSingleBlockMesh(c2d, cs, d, mpiRank);
    cs->msh->solveForJacobians();

    ///////////////////////////////////////////
    //Initialize Execution Loop and RK Method//
    ///////////////////////////////////////////
    AbstractRK *rk;
    rk = new TVDRK3(cs);

    IF_RANK0{
	cout << cs->msh->x_max[0] << endl;
	cout << cs->msh->x_max[1] << endl;
	cout << cs->msh->x_max[2] << endl;
    }

    ///////////////////////////////
    //Set flow initial conditions//
    ///////////////////////////////

    bool fromRestart = true;

    if(!fromRestart){
        FOR_Z_YPEN{
            FOR_Y_YPEN{
                FOR_X_YPEN{

                    int ip = GETMAJIND_YPEN;

		    int ii = GETGLOBALXIND_YPEN;		
		    int jj = GETGLOBALYIND_YPEN;		
		    int kk = GETGLOBALZIND_YPEN;

		    double x  = cs->msh->x[ip];
		    double x0 = 2.0*cs->msh->x_max[0]/8.0; 		
		    double y  = cs->msh->y[ip]; 		
		    double y0 = 1.2*cs->msh->x_max[1]/2.0; 		
		    double z  = cs->msh->z[ip]; 		
		    double z0 = cs->msh->x_max[2]/2.0; 		

		    double r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0);

                    cs->rho0[ip] = 1.0;//0.4 + 3.0*(x/20.0);
                    cs->p0[ip]   = (1.0 + 2.0*exp(-r2/0.001))/cs->ig->gamma;
                    cs->U0[ip]   = 0.2;
                    cs->V0[ip]   = 0.0;
                    cs->W0[ip]   = 0.0;
                }
            }
        }
    }else{

	string filename = "SolutionDump.15000";
	int timestep_start = 15000;

	cs->timeStep = timestep_start;

	MPI_File fh;
	MPI_Offset disp, filesize;

	MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	disp = 0;

	cs->c2d->readVar(fh, disp, 1, cs->rho0);
	cs->c2d->readVar(fh, disp, 1, cs->rhoU1);
	cs->c2d->readVar(fh, disp, 1, cs->rhoV1);
	cs->c2d->readVar(fh, disp, 1, cs->rhoW1);
	cs->c2d->readVar(fh, disp, 1, cs->rhoE1);

	MPI_File_close(&fh);

	FOR_XYZ_YPEN{
	    cs->U0[ip] = cs->ig->solveU(cs->rho0[ip], cs->rhoU1[ip]);
	    cs->V0[ip] = cs->ig->solveU(cs->rho0[ip], cs->rhoV1[ip]);
	    cs->W0[ip] = cs->ig->solveU(cs->rho0[ip], cs->rhoW1[ip]);
	    cs->p0[ip] = cs->ig->solvep(cs->rho0[ip], cs->rhoE1[ip], cs->U0[ip], cs->V0[ip], cs->W0[ip]);
	}

    }


    rk->executeSolverLoop();  

    //Now lets kill MPI
    MPI_Finalize();

    return 0;
}









