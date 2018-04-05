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
#include "Domain.hpp"
#include "BC.hpp" 
#include "Derivatives.hpp"

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
 
    /////////////////////////
    //Initialize the Domain//
    /////////////////////////
    int    Nx = 64,
           Ny = 64,
           Nz = 64;
    double Lx = 172.0,
           Ly = 129.0,
           Lz = 86.0;;
    Domain *d = new Domain(Nx, Ny, Nz, Lx, Ly, Lz, mpiRank);

    ///////////////////////////
    //Boundary Condition Info//
    ///////////////////////////
    BC::BCType bcXType = BC::PERIODIC_SOLVE;
    BC::BCType bcYType = BC::DIRICHLET_SOLVE;
    BC::BCType bcZType = BC::PERIODIC_SOLVE;

    BC::BCKind bcX0 = BC::PERIODIC;
    BC::BCKind bcX1 = BC::PERIODIC;
    BC::BCKind bcY0 = BC::SPONGE;
    BC::BCKind bcY1 = BC::SPONGE;
    BC::BCKind bcZ0 = BC::PERIODIC;
    BC::BCKind bcZ1 = BC::PERIODIC;

    bool periodicBC[3];
    BC *bc = new BC(bcXType, bcX0, bcX1,
                    bcYType, bcY0, bcY1,
                    bcZType, bcZ0, bcZ1,
		    periodicBC, mpiRank);



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
    d->setPencilDecompInfo(c2d->xSize,  c2d->ySize,  c2d->zSize,
				  c2d->xStart, c2d->yStart, c2d->zStart,
				  c2d->xEnd,   c2d->yEnd,   c2d->zEnd);
     

    //Initialize derivative objects for testing...
    Derivatives *derivX = new Derivatives(d, bc->bcXType, Derivatives::DIRX);
    Derivatives *derivY = new Derivatives(d, bc->bcYType, Derivatives::DIRY);
    Derivatives *derivZ = new Derivatives(d, bc->bcZType, Derivatives::DIRZ);


    double *u1, *u2, *u3;
    c2d->allocX(u1);
    c2d->allocX(u2);
    c2d->allocX(u3);

    c2d->deallocXYZ(u1);
    c2d->deallocXYZ(u2);
    c2d->deallocXYZ(u3);

    //Now lets kill MPI
    MPI_Finalize();

    return 0;
}









