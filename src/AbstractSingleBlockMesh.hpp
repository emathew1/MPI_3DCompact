#ifndef _CABSTRACTSINGLEBLOCKMESHH_
#define _CABSTRACTSINGLEBLOCKMESHH_

#include <iostream>
#include <fstream>

#include "AbstractCSolver.hpp"
#include "C2Decomp.hpp"
#include "Domain.hpp"
#include "Derivatives.hpp"
#include "Adt.hpp"

//There's some cyclic dependency going on here...
class AbstractCSolver;

class AbstractSingleBlockMesh{

    public:

	C2Decomp *c2d;
	AbstractCSolver *cs;
	Domain *d;

	int mpiRank;
	double *x, *y, *z;

	double *J;
	double *J11, *J12, *J13;
	double *J21, *J22, *J23;
	double *J31, *J32, *J33;

	double x_min[3], x_max[3];
	double max_xi, max_eta, max_zta;
	int Nx, Ny, Nz;

        int pxSize[3], pySize[3], pzSize[3];
        int pxStart[3], pyStart[3], pzStart[3];
        int pxEnd[3], pyEnd[3], pzEnd[3];

        double periodicXTranslation[3];
        double periodicYTranslation[3];
        double periodicZTranslation[3];

	bool periodicBCX, periodicBCY, periodicBCZ;
	bool transPeriodicX, transPeriodicY, transPeriodicZ;
	bool interPeriodicX, interPeriodicY, interPeriodicZ;

	Derivatives *derivX, *derivY, *derivZ;

	Adt<double> *adt;

	void allocateForMesh();
	void solveForJacobians();
	void dumpGrid();
	void getOrderedBlockCoordinates(int ip, int jp, int kp, double *x_halo, double *y_halo, double *z_halo, double box_p[8][3]);
	void getOrderedBlockXiCoordinates(int ip, int jp, int kp, double box_pxi[8][3]);
	void generateCoordinateHaloArrays(double *&x_halo, double *&y_halo, double *&z_halo);
	int findCVForPoint(double p[3], double *x_halo, double *y_halo, double *z_halo);
	void initMeshADT();
	void handlePeriodicStuff();
	

	//This is the function that makes this class abstract, mesh is never read or generated within this class...
	virtual void getMesh() = 0;

};


#endif
