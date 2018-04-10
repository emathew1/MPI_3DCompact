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
    int    Nx = 20,
           Ny = 20,
           Nz = 20;
    double Lx = 19.0,
           Ly = 19.0,
           Lz = 19.0;;
    Domain *d = new Domain(Nx, Ny, Nz, Lx, Ly, Lz, mpiRank);

    ///////////////////////////
    //Boundary Condition Info//
    ///////////////////////////
    BC::BCType bcXType = BC::DIRICHLET_SOLVE;
    BC::BCType bcYType = BC::DIRICHLET_SOLVE;
    BC::BCType bcZType = BC::DIRICHLET_SOLVE;

    BC::BCKind bcX0 = BC::SPONGE;
    BC::BCKind bcX1 = BC::SPONGE;
    BC::BCKind bcY0 = BC::SPONGE;
    BC::BCKind bcY1 = BC::SPONGE;
    BC::BCKind bcZ0 = BC::SPONGE;
    BC::BCKind bcZ1 = BC::SPONGE;

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
    d->setPencilDecompInfo(c2d);

    //This is info needed to use the macros...
    int pxSize[3], pySize[3], pzSize[3];
    int pxStart[3], pyStart[3], pzStart[3];
    int pxEnd[3], pyEnd[3], pzEnd[3];
    d->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);
     

    //Initialize derivative objects for testing...
    Derivatives *derivX = new Derivatives(d, bc->bcXType, Derivatives::DIRX);
    Derivatives *derivY = new Derivatives(d, bc->bcYType, Derivatives::DIRY);
    Derivatives *derivZ = new Derivatives(d, bc->bcZType, Derivatives::DIRZ);


    double *u1, *v1, *w1;
    double *du1dx;
    double *v2, *dv2dy, *dv1dy;
    double *w3, *dw3dz, *dw1dz;

    c2d->allocX(u1); 
    c2d->allocX(du1dx);

    c2d->allocX(v1);
    c2d->allocX(dv1dy);
    
    c2d->allocX(w1);
    c2d->allocX(dw1dz);

    c2d->allocY(v2);
    c2d->allocY(dv2dy);

    c2d->allocZ(w3);
    c2d->allocZ(dw3dz);
    

    FOR_X_XPEN{
	FOR_Y_XPEN{
	    FOR_Z_XPEN{

		int ip = GETMAJIND_XPEN;
		int ii = GETGLOBALXIND_XPEN; 
		int jj = GETGLOBALYIND_XPEN; 
		int kk = GETGLOBALZIND_XPEN; 

		u1[ip] = (double)ii;
		v1[ip] = (double)jj;
		w1[ip] = (double)kk;

	    }
	}
    }

    //Getting the X derivative
    derivX->calc1stDerivField(u1, du1dx); 
 
    //Getting the Y derivative
    c2d->transposeX2Y_MajorIndex(v1, v2);
    derivY->calc1stDerivField(v2, dv2dy);
    c2d->transposeY2X_MajorIndex(dv2dy, dv1dy);

    //Getting the Z derivative
    c2d->transposeX2Y_MajorIndex(w1, v2);
    c2d->transposeY2Z_MajorIndex(v2, w3);
    derivZ->calc1stDerivField(w3, dw3dz); 
    c2d->transposeZ2Y_MajorIndex(dw3dz, v2);
    c2d->transposeY2Z_MajorIndex(v2, dw1dz);

    IF_RANK0{
        FOR_X_XPEN{
	    FOR_Y_XPEN{
	        FOR_Z_XPEN{
		    int ip = GETMAJIND_XPEN;
		    cout << dv1dy[ip] << " "; 
	        }
		cout << endl;
	    }
	    cout << endl;
        }
	cout << endl;
    }


    c2d->deallocXYZ(u1);
    c2d->deallocXYZ(v1);
    c2d->deallocXYZ(v2);
    c2d->deallocXYZ(w1);
    c2d->deallocXYZ(w3);

    c2d->deallocXYZ(du1dx);
    c2d->deallocXYZ(dv1dy);
    c2d->deallocXYZ(dw1dz);
    c2d->deallocXYZ(dw3dz);
    c2d->deallocXYZ(dv2dy);
    

    //Now lets kill MPI
    MPI_Finalize();

    return 0;
}









