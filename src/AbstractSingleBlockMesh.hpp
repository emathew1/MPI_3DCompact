#ifndef _CABSTRACTSINGLEBLOCKMESHH_
#define _CABSTRACTSINGLEBLOCKMESHH_

#include "Derivatives.hpp"
#include "Adt.hpp"

class AbstractSingleBlockMesh{

    public:

	int mpiRank;
	double *x, *y, *z;

	double *J;
	double *J11, *J12, *J13;
	double *J21, *J22, *J23;
	double *J31, *J32, *J33;

	double max_xi, max_eta, max_zta;
	int Nx, Ny, Nz;

        int pxSize[3], pySize[3], pzSize[3];
        int pxStart[3], pyStart[3], pzStart[3];
        int pxEnd[3], pyEnd[3], pzEnd[3];

        double periodicXTranslation[3];
        double periodicYTranslation[3];
        double periodicZTranslation[3];

	bool periodicX, periodicY, periodicZ;

	Derivatives *derivX, *derivY, *derivZ;

	Adt<double> *adt;

	//Each RK Class needs to have these functions to overwrite the pure virtual ones
	virtual void solveForJacobians() = 0;

	virtual void generateCoordinateHaloArrays(double *&x_halo, double *&y_halo, double *&z_halo) = 0;
	
	virtual void getOrderedBlockCoordinates(int ip, int jp, int kp, double *x_halo, double *y_halo, double *z_halo, double box_p[8][3]) = 0;
	virtual void getOrderedBlockXiCoordinates(int ip, int jp, int kp, double box_pxi[8][3]) = 0;

	virtual int findCVForPoint(double p[3], double *x_halo, double *y_halo, double *z_halo) = 0;

};

#endif
