#ifndef _CURVILINEARINTEROLATORH_
#define _CURVILINEARINTEROLATORH_

#include "Macros.hpp"
#include "Utils.hpp"
#include "Domain.hpp"
#include "BC.hpp"
#include "AbstractSingleBlockMesh.hpp"
#include <vector>

class CurvilinearInterpolator{

  public:

    int pxSize[3], pySize[3], pzSize[3];
    int pxStart[3], pyStart[3], pzStart[3];
    int pxEnd[3], pyEnd[3], pzEnd[3];

    AbstractCSolver *cs;
    double (*pointList)[3];
    vector<int> icvList;
    int localPointFoundCount, globalPointFoundCount;

    //Weights for the components of the box to interpolate
    double (*Ni)[8];

    CurvilinearInterpolator(AbstractCSolver *cs, double (*pointList)[3], int Npoints){
	this->cs = cs;;
	this->pointList = pointList;

	cs->dom->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);

	//Start getting the points to interpolate
	double *xHalo, *yHalo, *zHalo;
	cs->msh->generateCoordinateHaloArrays(xHalo, yHalo, zHalo);
	
	localPointFoundCount = 0;
	for(int ip = 0; ip < Npoints;  ip++){

	    double p[3] = {pointList[ip][0], pointList[ip][1], pointList[ip][2]};
	
	    int icv = cs->msh->findCVForPoint(p, xHalo, yHalo, zHalo);
		
  	    if(icv != -1){
		localPointFoundCount++;
		icvList.push_back(icv);
	    }
	}

	MPI_Allreduce(&localPointFoundCount, &globalPointFoundCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 

        if(cs->mpiRank == 0){
	    cout << " > Found " << globalPointFoundCount << " of " << Npoints << " points for interpolation in the domain." << endl;
	}

	//allocate the local weights
	Ni = new double[localPointFoundCount++][8];

	//For each of local found points
	for(int ii = 0; ii < localPointFoundCount; ii++){

	    double f[8][3];
		
	    //To confuse things even more, box is oriented differently
	    //point 1 ->  1  1  1, box_p index -> 7
	    //point 2 -> -1  1  1, box_p index -> 3
	    //point 3 -> -1  1 -1, box_p index -> 2
	    //point 4 ->  1  1 -1, box_p index -> 6
	    //point 5 ->  1 -1  1, box_p index -> 5
	    //point 6 -> -1 -1  1, box_p index -> 1
	    //point 7 -> -1 -1 -1, box_p index -> 0
	    //point 8 ->  1 -1 -1, box_p index -> 4

	    //back out the ip, jp, kp based off index
	    int icv = icvList[ii];

	    //y-pencil major index oriented
	    int jp =  icv%pySize[1];
	    int kp = (icv/pySize[1])%pySize[2];
	    int ip =  icv/(pySize[2]*pySize[1]);

	    double box_p[8][3];

	    //Get the coordinates for this point
	    cs->msh->getOrderedBlockCoordinates(ip, jp, kp, xHalo, yHalo, zHalo, box_p); 

	    //reorder to make things easier to program
	    double xp[8][3];
	    FOR_I3{
		xp[0][i] = box_p[7][i];
		xp[1][i] = box_p[3][i];
		xp[2][i] = box_p[2][i];
		xp[3][i] = box_p[6][i];
		xp[4][i] = box_p[5][i];
		xp[5][i] = box_p[1][i];
		xp[6][i] = box_p[0][i];
		xp[7][i] = box_p[4][i];
	    }

	    //Get the f values for x, y, & z
	    FOR_I3{
		f[0][i] =  xp[0][i] + xp[1][i] + xp[2][i] + xp[3][i] + xp[4][i] + xp[5][i] + xp[6][i] + xp[7][i];
		f[1][i] =  xp[0][i] - xp[1][i] - xp[2][i] + xp[3][i] + xp[4][i] - xp[5][i] - xp[6][i] + xp[7][i];
		f[2][i] =  xp[0][i] + xp[1][i] + xp[2][i] + xp[3][i] - xp[4][i] - xp[5][i] - xp[6][i] - xp[7][i];
		f[3][i] =  xp[0][i] + xp[1][i] - xp[2][i] - xp[3][i] + xp[4][i] + xp[5][i] - xp[6][i] - xp[7][i];
		f[4][i] =  xp[0][i] - xp[1][i] - xp[2][i] + xp[3][i] - xp[4][i] + xp[5][i] + xp[6][i] - xp[7][i];
		f[5][i] =  xp[0][i] + xp[1][i] - xp[2][i] - xp[3][i] - xp[4][i] - xp[5][i] + xp[6][i] + xp[7][i];
		f[6][i] =  xp[0][i] - xp[1][i] + xp[2][i] - xp[3][i] + xp[4][i] - xp[5][i] + xp[6][i] - xp[7][i];
		f[7][i] =  xp[0][i] - xp[1][i] + xp[2][i] - xp[3][i] - xp[4][i] + xp[5][i] - xp[6][i] + xp[7][i];
	    }

	    for(int i = 0; i < 8; i++){
		for(int j = 0; j < 3; j++){
		    f[i][j] /= (1.0/8.0);
		}
	    }

	
	}

    }


};



#endif
