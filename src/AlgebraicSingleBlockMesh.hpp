#ifndef _CALGEBRAICSINGLEBLOCKMESHH_
#define _CALGEBRAICSINGLEBLOCKMESHH_

#include "Utils.hpp"
#include "Derivatives.hpp"
#include "AbstractSingleBlockMesh.hpp"

class AlgebraicSingleBlockMesh:public AbstractSingleBlockMesh{

    public:

	double a1, a2, a3;
	double b2, b3;
	double c1, c2, c3;

        int pxSize[3], pySize[3], pzSize[3]; 
        int pxStart[3], pyStart[3], pzStart[3];
        int pxEnd[3], pyEnd[3], pzEnd[3];

	AlgebraicSingleBlockMesh(){
	}


        void solveForJacobians();

};


void AlgebraicSingleBlockMesh::solveForJacobians(){
 
}

#endif
