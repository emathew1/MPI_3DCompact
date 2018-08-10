#ifndef _CURVILINEARINTEROLATORH_
#define _CURVILINEARINTEROLATORH_

#include "Macros.hpp"
#include "Utils.hpp"
#include "Domain.hpp"
#include "BC.hpp"
#include "AbstractSingleBlockMesh.hpp"
#include "Adt.hpp"

class CurvilinearInterpolator{

  public:

    Domain *c2d;
    BC *bc;
    AbstractSingleBlockMesh *msh;
    double (*pointList)[3];
    Adt<double> adt;


    CurvilinearInterpolator(Domain *c2d, BC *bc, AbstractSingleBlockMesh *msh, double (*pointList)[3]){
	this->c2d = c2d;
	this->BC  = bc;
	this->msh = msh;
	this->pointList = pointList;
    }


};



#endif
