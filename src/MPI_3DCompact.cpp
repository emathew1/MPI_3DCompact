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
    double CFL       = 0.8;
    int maxTimeStep  = 250;
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
    BC::BCType bcXType = BC::PERIODIC_SOLVE;
    BC::BCType bcYType = BC::PERIODIC_SOLVE;
    BC::BCType bcZType = BC::PERIODIC_SOLVE;

    BC::BCKind bcX0 = BC::PERIODIC;
    BC::BCKind bcX1 = BC::PERIODIC;
    BC::BCKind bcY0 = BC::PERIODIC;
    BC::BCKind bcY1 = BC::PERIODIC;
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
    int    Nx = 100,
           Ny = 100,
           Nz = 100;
    double Lx = 1.0,
           Ly = 1.0,
           Lz = 1.0;;
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
    double alphaF  = 0.375;
    double mu_ref  = 0.00375;
    bool useTiming = false;
    AbstractCSolver *cs;
    cs = new CurvilinearCSolver(c2d, d, bc, ts, alphaF, mu_ref, useTiming);
    
    AbstractSingleBlockMesh *msh;
    msh = new AlgebraicSingleBlockMesh(c2d, cs, d, mpiRank);
    msh->solveForJacobians();


    getRange(msh->J, "J", pySize[0], pySize[1], pySize[2], mpiRank);
 
    double *J11;
    c2d->allocY(J11);
    FOR_XYZ_YPEN{
	J11[ip] = msh->J[ip]*msh->J11[ip];
    }

    getRange(msh->J11, "J11", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(msh->J12, "J12", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(msh->J13, "J13", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(msh->J21, "J21", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(msh->J22, "J22", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(msh->J23, "J23", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(msh->J31, "J31", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(msh->J32, "J32", pySize[0], pySize[1], pySize[2], mpiRank);
    getRange(msh->J33, "J33", pySize[0], pySize[1], pySize[2], mpiRank);

    ///////////////////////////////////////////
    //Initialize Execution Loop and RK Method//
    ///////////////////////////////////////////
    AbstractRK *rk;
    rk = new TVDRK3(cs);


    ///////////////////////////////
    //Set flow initial conditions//
    ///////////////////////////////
    FOR_Z_YPEN{
        FOR_Y_YPEN{
            FOR_X_YPEN{

                int ip = GETMAJIND_YPEN;

		int ii = GETGLOBALXIND_YPEN;		
		int jj = GETGLOBALYIND_YPEN;		
		int kk = GETGLOBALZIND_YPEN;

		double x  = d->x[ii];
		double x0 = Lx/2; 		
		double y  = d->y[jj]; 		
		double y0 = Ly/2; 		
		double z  = d->z[kk]; 		
		double z0 = Lz/2; 		

		double r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0);

                cs->rho0[ip] = 1.0;
                cs->p0[ip]   = (1.0 + 10.0*exp(-r2/0.001))/cs->ig->gamma;
                cs->U0[ip]   = 0.0;
                cs->V0[ip]   = 0.0;
                cs->W0[ip]   = 0.0;
            }
        }
    }
   
//    rk->executeSolverLoop();  

    //Now lets kill MPI
    MPI_Finalize();

    return 0;
}









