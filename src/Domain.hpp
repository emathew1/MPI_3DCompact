#ifndef _DOMAINH_
#define _DOMAINH_

#include <iostream>
#include "Macros.hpp"
#include "C2Decomp.hpp"

class Domain{

    public:

	//Global variables are denoted with a g
        int gNx, gNy, gNz, gN;
        double gLx, gLy, gLz;

	//Local varibles...
	int pxSize[3], pySize[3], pzSize[3]; //pencil x, y, z, Sizes in 3 dimensions
	int pxStart[3], pyStart[3], pzStart[3];
	int pxEnd[3], pyEnd[3], pzEnd[3];

	//These are going to be based on the global indices
        double *x, *y, *z;

	double dx, dy, dz;

	int mpiRank;


    Domain(int Nx, int Ny, int Nz, double Lx, double Ly, double Lz, int mpiRank){
	
	gNx = Nx;
	gNy = Ny;
	gNz = Nz;
	gN = gNx*gNy*gNz;

	gLx = Lx;
	gLy = Ly;
	gLz = Lz;

	x = new double[gNx];	
	y = new double[gNy];	
	z = new double[gNz];

        for(int ip = 0; ip < gNx; ip++){
	    x[ip] = (((double)ip)/((double)gNx - 1.0))*gLx;
	}	

        for(int jp = 0; jp < gNy; jp++){
	    y[jp] = (((double)jp)/((double)Ny - 1.0))*gLy;
	}	

        for(int kp = 0; kp < gNz; kp++){
	    z[kp] = (((double)kp)/((double)Nz - 1.0))*gLz;
	}	

	dx = x[1]-x[0];
	dy = y[1]-y[0];
	dz = z[1]-z[0];

	if(!mpiRank){
	    std::cout << " > Global Domain initialization..." << std::endl;
	    std::cout << " > Domain: " << gLx << "x" << gLy << "x" << gLz << std::endl;
	    std::cout << " > Mesh: " << gNx << "x" << gNy << "x" << gNz << std::endl;
	    std::cout << " > Total Points: " << gNx*gNy*gNz << std::endl; 
	}
    }
  
    void setPencilDecompInfo(C2Decomp *c2d){
	FOR_I3 pxSize[i] = c2d->xSize[i];
	FOR_I3 pySize[i] = c2d->ySize[i];
	FOR_I3 pzSize[i] = c2d->zSize[i];

	FOR_I3 pxStart[i] = c2d->xStart[i];
	FOR_I3 pyStart[i] = c2d->yStart[i];
	FOR_I3 pzStart[i] = c2d->zStart[i];

	FOR_I3 pxEnd[i] = c2d->xEnd[i];
	FOR_I3 pyEnd[i] = c2d->yEnd[i];
	FOR_I3 pzEnd[i] = c2d->zEnd[i];
	
    } 

    void getPencilDecompInfo(int xSize[3], int ySize[3], int zSize[3], int xStart[3], int yStart[3], int zStart[3], int xEnd[3], int yEnd[3], int zEnd[3]){

	FOR_I3 xSize[i] = pxSize[i];
	FOR_I3 ySize[i] = pySize[i];
	FOR_I3 zSize[i] = pzSize[i];

	FOR_I3 xStart[i] = pxStart[i];
	FOR_I3 yStart[i] = pyStart[i];
	FOR_I3 zStart[i] = pzStart[i];

	FOR_I3 xEnd[i] = pxEnd[i];
	FOR_I3 yEnd[i] = pyEnd[i];
	FOR_I3 zEnd[i] = pzEnd[i];
	
    }

};



#endif
