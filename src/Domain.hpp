#ifndef _DOMAINH_
#define _DOMAINH_

#include <iostream>

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

};

#endif
