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
        cout << "-------------------" << endl;
    	cout << " C2Decomp Testing " << endl;
    	cout << "-------------------" << endl;
    	cout << endl;

    }
 
    int nx = 20, ny = 15, nz = 10;
    int pRow = 0, pCol = 0;
    bool periodicBC[3] = {true, true, true};

    if(!mpiRank) cout << "initializing " << endl;
    C2Decomp *c2d;
    c2d = new C2Decomp(nx, ny, nz, pRow, pCol, periodicBC);
    if(!mpiRank) cout << "done initializing " << endl;
    
    bool errorFlag, errorFlagGlobal;

    int m = 1;
    double data1[nz][ny][nx];
    for(int kp = 0; kp < nz; kp++){
	for(int jp = 0; jp < ny; jp++){
	    for(int ip = 0; ip < nx; ip++){
		data1[kp][jp][ip] = (double)m;
		m++;
	    }
	}
    }

    double xSize[3], ySize[3], zSize[3];
    xSize[0] = c2d->xSize[0];
    xSize[1] = c2d->xSize[1];
    xSize[2] = c2d->xSize[2];
    ySize[0] = c2d->ySize[0];
    ySize[1] = c2d->ySize[1];
    ySize[2] = c2d->ySize[2];
    zSize[0] = c2d->zSize[0];
    zSize[1] = c2d->zSize[1];
    zSize[2] = c2d->zSize[2];

    double *u1, *u2, *u3;


    c2d->allocX(u1);
    c2d->allocY(u2);
    c2d->allocZ(u3);

    for(int kp = 0; kp < xSize[2]; kp++){
	for(int jp = 0; jp < xSize[1]; jp++){
	    for(int ip = 0; ip < xSize[0]; ip++){
		int ii = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip;
		u1[ii] = data1[c2d->xStart[2]+kp][c2d->xStart[1]+jp][c2d->xStart[0]+ip];
	    }
	}
    }

    double t1, t2, t3;
    t1 = MPI_Wtime();
    c2d->transposeX2Y_MajorIndex(u1, u2);
    t2 = MPI_Wtime();

    if(mpiRank == 0){
	printf( "X2Y Elapsed time is %f\n", t2 - t1 );
    }

    //Testing transposition

    errorFlag = false;
    errorFlagGlobal = false;
    for(int kp = 0; kp < ySize[2]; kp++){
	for(int jp = 0; jp < ySize[1]; jp++){
	    for(int ip = 0; ip < ySize[0]; ip++){
	  	int ii = ip*ySize[2]*ySize[1] + kp*ySize[1] + jp;
		double temp = u2[ii];
		double temp1= data1[c2d->yStart[2]+kp][c2d->yStart[1]+jp][c2d->yStart[0]+ip];  
		if(fabs(temp-temp1)>1E-16){
		    cout << "Error in blocking X2Y tranposition" << endl;
		}
	    }
	}
    }


    t1 = MPI_Wtime();
    c2d->transposeY2Z_MajorIndex(u2, u3);
    t2 = MPI_Wtime();
    if(mpiRank == 0){
	printf( "Y2Z Elapsed time is %f\n", t2 - t1 );
    }


    //Testing transposition
    for(int kp = 0; kp < zSize[2]; kp++){
	for(int jp = 0; jp < zSize[1]; jp++){
	    for(int ip = 0; ip < zSize[0]; ip++){
		int ii = jp*zSize[2]*zSize[0] + ip*zSize[2] + kp;
		double temp = u3[ii];
		double temp1= data1[c2d->zStart[2]+kp][c2d->zStart[1]+jp][c2d->zStart[0]+ip]; 
		if(fabs(temp-temp1)>1E-16){
		    cout << "Error in blocking Y2Z tranposition" << endl;
		}
 	    }
  	}
    }
    

    t1 = MPI_Wtime();
    c2d->transposeZ2Y_MajorIndex(u3, u2);
    t2 = MPI_Wtime();
    if(mpiRank == 0){
	printf( "Z2Y Elapsed time is %f\n", t2 - t1 );
    }


    //Testing transposition
    for(int kp = 0; kp < ySize[2]; kp++){
	for(int jp = 0; jp < ySize[1]; jp++){
	    for(int ip = 0; ip < ySize[0]; ip++){
	  	int ii = ip*ySize[2]*ySize[1] + kp*ySize[1] + jp;
		double temp = u2[ii];
		double temp1= data1[c2d->yStart[2]+kp][c2d->yStart[1]+jp][c2d->yStart[0]+ip]; 
		if(fabs(temp-temp1)>1E-16){
		    cout << "Error in blocking Z2Y tranposition" << endl;
		}
	    }
 	}
    }


    t1 = MPI_Wtime();
    c2d->transposeY2X_MajorIndex(u2, u1);
    t2 = MPI_Wtime();
    if(mpiRank == 0){
	printf( "Y2X Elapsed time is %f\n", t2 - t1 );
    }

    if(mpiRank==0){
      //Testing transposition
      for(int kp = 0; kp < xSize[2]; kp++){
	  for(int jp = 0; jp < xSize[1]; jp++){
	      for(int ip = 0; ip < xSize[0]; ip++){
	  	  int ii = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip;
		  double temp = u1[ii];
		  double temp1= data1[c2d->xStart[2]+kp][c2d->xStart[1]+jp][c2d->xStart[0]+ip]; 
		  if(fabs(temp-temp1)>1E-16){
		      cout << "Error in blocking Y2X transposition" << endl;
		  }
	      }
 	  }
       }
    }

    //allocate new buffers for non-blocking comms
    double *sbuf = new double[c2d->decompBufSize];
    double *rbuf = new double[c2d->decompBufSize];

    MPI_Request x2yHandle;

    t1 = MPI_Wtime();
    c2d->transposeX2Y_MajorIndex_Start(x2yHandle, u1, u2, sbuf, rbuf);
    t2 = MPI_Wtime();
    c2d->transposeX2Y_MajorIndex_Wait(x2yHandle, u1, u2, sbuf, rbuf);
    t3 = MPI_Wtime();
    if(mpiRank == 0){
	printf( "X2Y Nonblocking Start Elapsed time is %f\n", t2 - t1 );
	printf( "X2Y Nonblocking Wait Elapsed time is %f\n", t3 - t2 );
    }


    //Testing transposition
    for(int kp = 0; kp < ySize[2]; kp++){
	for(int jp = 0; jp < ySize[1]; jp++){
	    for(int ip = 0; ip < ySize[0]; ip++){
	  	int ii = ip*ySize[1]*ySize[2] + kp*ySize[1] + jp;
		double temp = u2[ii];
		double temp1= data1[c2d->yStart[2]+kp][c2d->yStart[1]+jp][c2d->yStart[0]+ip]; 
		if(fabs(temp-temp1)>1E-16){
		    cout << "Error in nonblocking X2Y transposition" << endl;
		}
	    }
        }
    }

    MPI_Request y2zHandle;
    t1 = MPI_Wtime();
    c2d->transposeY2Z_MajorIndex_Start(y2zHandle, u2, u3, sbuf, rbuf);
    t2 = MPI_Wtime();
    c2d->transposeY2Z_MajorIndex_Wait(y2zHandle, u2, u3, sbuf, rbuf);
    t3 = MPI_Wtime();
    if(mpiRank == 0){
	printf( "Y2Z Nonblocking Start Elapsed time is %f\n", t2 - t1 );
	printf( "Y2Z Nonblocking Wait Elapsed time is %f\n", t3 - t2 );
    }


    //Testing transposition
    for(int kp = 0; kp < zSize[2]; kp++){
        for(int jp = 0; jp < zSize[1]; jp++){
	    for(int ip = 0; ip < zSize[0]; ip++){
		int ii = jp*zSize[2]*zSize[0] + ip*zSize[2] + kp;
		double temp = u3[ii];
		double temp1= data1[c2d->zStart[2]+kp][c2d->zStart[1]+jp][c2d->zStart[0]+ip];
		if(fabs(temp-temp1)>1E-16){  
		    cout << "Error in nonblocking Y2Z transposition" << endl;
		}
 	    }
  	}
    }


    MPI_Request z2yHandle;
    t1 = MPI_Wtime();
    c2d->transposeZ2Y_MajorIndex_Start(z2yHandle, u3, u2, sbuf, rbuf);
    t2 = MPI_Wtime();
    c2d->transposeZ2Y_MajorIndex_Wait(z2yHandle, u3, u2, sbuf, rbuf);
    t3 = MPI_Wtime();
    if(mpiRank == 0){
	printf( "Z2Y Nonblocking Start Elapsed time is %f\n", t2 - t1 );
	printf( "Z2Y Nonblocking Wait Elapsed time is %f\n", t3 - t2 );
    }


    //Testing transposition
    for(int kp = 0; kp < ySize[2]; kp++){
	for(int jp = 0; jp < ySize[1]; jp++){
	    for(int ip = 0; ip < ySize[0]; ip++){
	  	int ii = ip*ySize[1]*ySize[2] + kp*ySize[1] + jp;
		double temp = u2[ii];
		double temp1= data1[c2d->yStart[2]+kp][c2d->yStart[1]+jp][c2d->yStart[0]+ip]; 
		if(fabs(temp-temp1)>1E-16){  
		    cout << "Error in nonblocking Z2Y transposition" << endl;
	   	}
	    }
 	}
    }


    MPI_Request y2xHandle;
    t1 = MPI_Wtime();
    c2d->transposeY2X_MajorIndex_Start(y2xHandle, u2, u1, sbuf, rbuf);
    t2 = MPI_Wtime();
    c2d->transposeY2X_MajorIndex_Wait(y2xHandle, u2, u1, sbuf, rbuf);
    t3 = MPI_Wtime();
    if(mpiRank == 0){
	printf( "Y2X Nonblocking Start Elapsed time is %f\n", t2 - t1 );
	printf( "Y2X Nonblocking Wait Elapsed time is %f\n", t3 - t2 );
    }


    //Testing transposition
    for(int kp = 0; kp < xSize[2]; kp++){
	for(int jp = 0; jp < xSize[1]; jp++){
	    for(int ip = 0; ip < xSize[0]; ip++){
	  	int ii = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip;
		double temp = u1[ii];
		double temp1= data1[c2d->xStart[2]+kp][c2d->xStart[1]+jp][c2d->xStart[0]+ip]; 
		if(fabs(temp-temp1)>1E-16){  
		    cout << "Error in nonblocking Y2X transposition" << endl;
	   	}
	    }
 	}
    }



    delete[] sbuf;
    delete[] rbuf;
    c2d->deallocXYZ(u1);
    c2d->deallocXYZ(u2);
    c2d->deallocXYZ(u3);

    

    //Now lets kill MPI
    MPI_Finalize();

    return 0;
}









