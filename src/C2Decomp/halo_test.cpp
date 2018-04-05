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


double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main(int argc, char *argv[]){

    int ierr, totRank, mpiRank;

    //Initialize MPI
    ierr = MPI_Init( &argc, &argv);

    //Get the number of processes
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &totRank);

    //Get the local rank
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    if(!mpiRank){
        cout << endl;
        cout << "-----------------------" << endl;
    	cout << " C2Decomp Halo Testing " << endl;
    	cout << "-----------------------" << endl;
    	cout << endl;

    }
 
    int nx = 100,  ny = 150, nz = 131;
    int pRow = 0, pCol = 0;
    bool periodicBC[3] = {true, true, true};

    if(!mpiRank) cout << "initializing " << endl;
    C2Decomp *c2d;
    c2d = new C2Decomp(nx, ny, nz, pRow, pCol, periodicBC);
    if(!mpiRank) cout << "done initializing " << endl;

    int xSize[3] = {c2d->xSize[0], c2d->xSize[1], c2d->xSize[2]};
    int ySize[3] = {c2d->ySize[0], c2d->ySize[1], c2d->ySize[2]};
    int zSize[3] = {c2d->zSize[0], c2d->zSize[1], c2d->zSize[2]};

    
    double *u1=NULL, *u2=NULL, *u3=NULL;
    double *v1=NULL, *v2=NULL, *v3=NULL;
    double *w1=NULL, *w2=NULL, *w3=NULL;

    double *wk2=NULL, *wk3=NULL;
    double *uh1=NULL, *vh1=NULL, *wh1=NULL;
    double *uh2=NULL, *vh2=NULL, *wh2=NULL;
    double *uh3=NULL, *vh3=NULL, *wh3=NULL;
    double *div1=NULL, *div2=NULL, *div3=NULL, *div4=NULL;

    c2d->allocX(u1);
    c2d->allocX(v1);
    c2d->allocX(w1);
    c2d->allocX(div1);
    c2d->allocX(div2);
    c2d->allocX(div3);
    c2d->allocX(div4);

    c2d->allocY(u2);
    c2d->allocY(v2);
    c2d->allocY(w2);
    c2d->allocY(wk2);

    c2d->allocZ(u3);
    c2d->allocZ(v3);
    c2d->allocZ(w3);
    c2d->allocZ(wk3);


 
    //Initialize the u, v, and w fields with random variables
    srand(1);

    if(!mpiRank) cout << "Initializing random field... " << endl;
    for(int kp = 0; kp < xSize[2]; kp++){
	for(int jp = 0; jp < xSize[1]; jp++){
	    for(int ip = 0; ip < xSize[0]; ip++){
		int ii = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip;
		u1[ii] = fRand(-1.0,1.0);
		v1[ii] = fRand(-1.0,1.0);
		w1[ii] = fRand(-1.0,1.0);
	    }
	}
    }

    //Initialize the first divergence test to zero

    for(int ip = 0; ip < xSize[2]*xSize[1]*xSize[0]; ip++){
	div1[ip] = 0.0;
    }

    if(!mpiRank) cout << "Calculating Divergence in xpencil... " << endl;
    //Calculate the divergence in the x divergence along the x pencil
     for(int kp = 0; kp < xSize[2]; kp++){
	for(int jp = 0; jp < xSize[1]; jp++){
	    for(int ip = 1; ip < xSize[0]-1; ip++){
		int ii   = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip;
		int iip1 = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip + 1;
		int iim1 = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip - 1;
		div1[ii] = u1[iip1]-u1[iim1];
	    }
	}
    }


    //Move the things we need to the y-pencil
    c2d->transposeX2Y(v1,v2);
    c2d->transposeX2Y(div1, wk2);

    if(!mpiRank) cout << "Calculating Divergence in ypencil... " << endl;
    for(int kp = 0; kp < ySize[2]; kp++){
	for(int jp = 1; jp < ySize[1]-1; jp++){
	    for(int ip = 0; ip < ySize[0]; ip++){
		int ii   = kp*ySize[1]*ySize[0] + jp*ySize[0]     + ip;
		int iip1 = kp*ySize[1]*ySize[0] + (jp+1)*ySize[0] + ip;
		int iim1 = kp*ySize[1]*ySize[0] + (jp-1)*ySize[0] + ip;
		wk2[ii] += v2[iip1]-v2[iim1];
	    }
	}
    }
    
    //Move the things we need to the z-pencil
    c2d->transposeX2Y(w1, w2);
    c2d->transposeY2Z(w2, w3);
    c2d->transposeY2Z(wk2, wk3);

    if(!mpiRank) cout << "Calculating Divergence in zpencil... " << endl;
    for(int kp = 1; kp < zSize[2]-1; kp++){
	for(int jp = 0; jp < zSize[1]; jp++){
	    for(int ip = 0; ip < zSize[0]; ip++){
		int ii   = kp*zSize[1]*zSize[0]     + jp*zSize[0] + ip;
		int iip1 = (kp+1)*zSize[1]*zSize[0] + jp*zSize[0] + ip;
		int iim1 = (kp-1)*zSize[1]*zSize[0] + jp*zSize[0] + ip;
		wk3[ii] += w3[iip1]-w3[iim1];
	    }
	}
    }

    //Move the divergence back to the x-pencil
    c2d->transposeZ2Y(wk3, wk2);
    c2d->transposeY2X(wk2, div1);

    //////////////////////////////////////
    //Now lets calculate it using halo's//
    //////////////////////////////////////


    if(!mpiRank) cout << "Calculating Divergence using halo's in xpencil" << endl;
    for(int kp = 0; kp < xSize[2]; kp++){
	for(int jp = 0; jp < xSize[1]; jp++){
	    for(int ip = 1; ip < xSize[0]-1; ip++){
		int ii   = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip;
		int iip1 = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip + 1;
		int iim1 = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip - 1;
		div2[ii] = u1[iip1]-u1[iim1];
	    }
	}
    }

    int numLevel = 1;
    int ipencil  = 0;
    c2d->updateHalo(v1, vh1, numLevel, ipencil); 

    for(int kp = 0; kp < xSize[2]; kp++){
	for(int jp = 0; jp < xSize[1]; jp++){
	    for(int ip = 1; ip < xSize[0]-1; ip++){
		int ii    = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip;

		//some macros or functions to retrieve this information for a given pencil would be nice...
		int iihp1 = (kp+numLevel)*(xSize[1]+2*numLevel)*xSize[0] + (jp+numLevel+1)*xSize[0] + ip;
		int iihm1 = (kp+numLevel)*(xSize[1]+2*numLevel)*xSize[0] + (jp+numLevel-1)*xSize[0] + ip;
		div2[ii] += vh1[iihp1]-vh1[iihm1];
	    }
	}
    }

    c2d->updateHalo(w1, wh1, numLevel, ipencil); 

    for(int kp = 0; kp < xSize[2]; kp++){
	for(int jp = 0; jp < xSize[1]; jp++){
	    for(int ip = 1; ip < xSize[0]-1; ip++){
		int ii    = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip;

		//some macros or functions to retrieve this information for a given pencil would be nice...
		int iihp1 = (kp+numLevel+1)*(xSize[1]+2*numLevel)*xSize[0] + (jp+numLevel)*xSize[0] + ip;
		int iihm1 = (kp+numLevel-1)*(xSize[1]+2*numLevel)*xSize[0] + (jp+numLevel)*xSize[0] + ip;
		div2[ii] += wh1[iihp1]-wh1[iihm1];
	    }
	}
    }

    bool errorFlag = false;
    bool errorFlagGlobal =false;
    for(int kp = 1; kp < xSize[2]-1; kp++){
	for(int jp = 1; jp < xSize[1]-1; jp++){
	    for(int ip = 1; ip < xSize[0]-1; ip++){
		int ii    = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip;
		if(fabs(div1[ii] - div2[ii]) > 1E-16){
		    errorFlag = true;
		}
	    }
	    
	}
    }

    MPI_Allreduce(&errorFlag, &errorFlagGlobal, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

    if(!mpiRank){
        if(errorFlagGlobal){
	    cout << "Error! Looks like a problem in x-pencil halo!" << endl;
        }else{
	    cout << "Halo support for x-pencil looks good!" << endl;
	}
    }


    MPI_Barrier(MPI_COMM_WORLD);

    c2d->deallocXYZ(vh1);
    c2d->deallocXYZ(wh1);

    if(!mpiRank) cout << "Calculating Divergence using halo's in ypencil" << endl;

    c2d->transposeX2Y(u1, u2);
    c2d->transposeX2Y(v1, v2);
    c2d->transposeX2Y(w1, w2);   
    

    numLevel = 1;
    ipencil  = 1;
    c2d->updateHalo(u2, uh2, numLevel, ipencil);
    MPI_Barrier(MPI_COMM_WORLD);

    for(int kp = 0; kp < ySize[2]; kp++){
	for(int jp = 1; jp < ySize[1]-1; jp++){
	    for(int ip = 0; ip < ySize[0]; ip++){
		int ii   = kp*ySize[1]*ySize[0] + jp*ySize[0]     + ip;

		int iihp1 = (kp+numLevel)*ySize[1]*(ySize[0]+2*numLevel) + jp*(ySize[0]+2*numLevel) + ip + numLevel + 1;
		int iihm1 = (kp+numLevel)*ySize[1]*(ySize[0]+2*numLevel) + jp*(ySize[0]+2*numLevel) + ip + numLevel - 1;
		wk2[ii] = uh2[iihp1]-uh2[iihm1];
	    }
	}
    }


    for(int kp = 0; kp < ySize[2]; kp++){
	for(int jp = 1; jp < ySize[1]-1; jp++){
	    for(int ip = 0; ip < ySize[0]; ip++){
		int ii   = kp*ySize[1]*ySize[0] + jp*ySize[0] + ip;
		int iip1 = kp*ySize[1]*ySize[0] + (jp+1)*ySize[0] + ip;
		int iim1 = kp*ySize[1]*ySize[0] + (jp-1)*ySize[0] + ip;
		wk2[ii] += v2[iip1]-v2[iim1];
	    }
	}
    }

    numLevel = 1;
    ipencil  = 1;
    c2d->updateHalo(w2, wh2, numLevel, ipencil);


    for(int kp = 0; kp < ySize[2]; kp++){
	for(int jp = 1; jp < ySize[1]-1; jp++){
	    for(int ip = 0; ip < ySize[0]; ip++){
		int ii   = kp*ySize[1]*ySize[0] + jp*ySize[0]     + ip;

		int iihp1 = (kp+numLevel+1)*ySize[1]*(ySize[0]+2*numLevel) + jp*(ySize[0]+2*numLevel) + ip + numLevel;
		int iihm1 = (kp+numLevel-1)*ySize[1]*(ySize[0]+2*numLevel) + jp*(ySize[0]+2*numLevel) + ip + numLevel;
		wk2[ii] += wh2[iihp1]-wh2[iihm1];
	    }
	}
    }

    c2d->transposeY2X(wk2, div3);


    errorFlag = false;
    errorFlagGlobal = false;
    for(int kp = 1; kp < xSize[2]-1; kp++){
	for(int jp = 1; jp < xSize[1]-1; jp++){
	    for(int ip = 1; ip < xSize[0]-1; ip++){
		int ii    = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip;
		if(fabs(div1[ii] - div3[ii]) > 1E-16){
		    errorFlag = true;
		}
	    }
	    
	}
    }

    MPI_Allreduce(&errorFlag, &errorFlagGlobal, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

    if(!mpiRank){
        if(errorFlagGlobal){
	    cout << "Error! Looks like a problem in y-pencil halo!" << endl;
        }else{
	    cout << "Halo support for y-pencil looks good!" << endl;
	}
    }


    MPI_Barrier(MPI_COMM_WORLD);

    c2d->deallocXYZ(uh2); 
    c2d->deallocXYZ(wh2); 

    if(!mpiRank) cout << "Calculating Divergence using halo's in zpencil" << endl;

    c2d->transposeY2Z(u2, u3); 
    c2d->transposeY2Z(v2, v3); 
    c2d->transposeY2Z(w2, w3);

    numLevel = 1;
    ipencil  = 2;
    c2d->updateHalo(u3, uh3, numLevel, ipencil);

    for(int kp = 1; kp < zSize[2]-1; kp++){
	for(int jp = 0; jp < zSize[1]; jp++){
	    for(int ip = 0; ip < zSize[0]; ip++){
		int ii   = kp*zSize[1]*zSize[0]     + jp*zSize[0] + ip;

		int iihp1 = (kp)*(zSize[1]+2*numLevel)*(zSize[0]+2*numLevel) + (jp+numLevel)*(zSize[0]+2*numLevel) + ip+numLevel+1;
		int iihm1 = (kp)*(zSize[1]+2*numLevel)*(zSize[0]+2*numLevel) + (jp+numLevel)*(zSize[0]+2*numLevel) + ip+numLevel-1;
		wk3[ii] = uh3[iihp1]-uh3[iihm1];
	    }
	}
    }

    numLevel = 1;
    ipencil  = 2;
    c2d->updateHalo(v3, vh3, numLevel, ipencil);

    for(int kp = 1; kp < zSize[2]-1; kp++){
	for(int jp = 0; jp < zSize[1]; jp++){
	    for(int ip = 0; ip < zSize[0]; ip++){
		int ii   = kp*zSize[1]*zSize[0]     + jp*zSize[0] + ip;

		int iihp1 = (kp)*(zSize[1]+2*numLevel)*(zSize[0]+2*numLevel) + (jp+numLevel+1)*(zSize[0]+2*numLevel) + ip+numLevel;
		int iihm1 = (kp)*(zSize[1]+2*numLevel)*(zSize[0]+2*numLevel) + (jp+numLevel-1)*(zSize[0]+2*numLevel) + ip+numLevel;
		wk3[ii] += vh3[iihp1]-vh3[iihm1];
	    }
	}
    }

    for(int kp = 1; kp < zSize[2]-1; kp++){
	for(int jp = 0; jp < zSize[1]; jp++){
	    for(int ip = 0; ip < zSize[0]; ip++){
		int ii   = kp*zSize[1]*zSize[0]     + jp*zSize[0] + ip;
		int iip1 = (kp+1)*zSize[1]*zSize[0] + jp*zSize[0] + ip;
		int iim1 = (kp-1)*zSize[1]*zSize[0] + jp*zSize[0] + ip;

		wk3[ii] += w3[iip1]-w3[iim1];
	    }
	}
    }

    c2d->transposeZ2Y(wk3, wk2);
    c2d->transposeY2X(wk2, div4);

    errorFlag = false;
    errorFlagGlobal = false;
    for(int kp = 1; kp < xSize[2]-1; kp++){
	for(int jp = 1; jp < xSize[1]-1; jp++){
	    for(int ip = 1; ip < xSize[0]-1; ip++){
		int ii    = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip;
		if(fabs(div1[ii] - div4[ii]) > 1E-16){
		    errorFlag = true;
		}
	    }
	    
	}
    }

    MPI_Allreduce(&errorFlag, &errorFlagGlobal, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

    if(!mpiRank){
        if(errorFlagGlobal){
	    cout << "Error! Looks like a problem in z-pencil halo!" << endl;
        }else{
	    cout << "Halo support for z-pencil looks good!" << endl;
	}
    }

    MPI_Barrier(MPI_COMM_WORLD);

    c2d->deallocXYZ(uh3);
    c2d->deallocXYZ(vh3);
 

    //Clean up our allocated memory...
    c2d->deallocXYZ(u1);
    c2d->deallocXYZ(v1);
    c2d->deallocXYZ(w1);
    c2d->deallocXYZ(u2);
    c2d->deallocXYZ(v2);
    c2d->deallocXYZ(w2);
    c2d->deallocXYZ(u3);
    c2d->deallocXYZ(v3);
    c2d->deallocXYZ(w3);
    c2d->deallocXYZ(wk2);
    c2d->deallocXYZ(wk3);
    c2d->deallocXYZ(div1);
    c2d->deallocXYZ(div2);
    c2d->deallocXYZ(div3);
    c2d->deallocXYZ(div4);
   

    //Now lets kill MPI
    MPI_Finalize();



    return 0;
}









