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

    AbstractCSolver *cs;
    double (*pointList)[3];
    vector<int> icvList;
    int localPointFoundCount, globalPointFoundCount;

    CurvilinearInterpolator(AbstractCSolver *cs, double (*pointList)[3], int Npoints){
	this->cs = cs;;
	this->pointList = pointList;

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

		

    }


};



#endif
