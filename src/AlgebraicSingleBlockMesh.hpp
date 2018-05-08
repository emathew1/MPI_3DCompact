#ifndef _CALGEBRAICSINGLEBLOCKMESHH_
#define _CALGEBRAICSINGLEBLOCKMESHH_

#include "Macros.hpp"
#include "Utils.hpp"
#include "Derivatives.hpp"
#include "AbstractSingleBlockMesh.hpp"

class AlgebraicSingleBlockMesh:public AbstractSingleBlockMesh{

    public:

	C2Decomp *c2d;
	Domain *dom;

        int pxSize[3], pySize[3], pzSize[3]; 
        int pxStart[3], pyStart[3], pzStart[3];
        int pxEnd[3], pyEnd[3], pzEnd[3];

	AlgebraicSingleBlockMesh(C2Decomp *c2d, Domain *dom, Derivatives *derivX, Derivatives *derivY, Derivatives *derivZ, int mpiRank){

	    this->mpiRank = mpiRank;

	    this->c2d = c2d;
	    this->dom = dom;
	    this->derivX = derivX;
	    this->derivY = derivY;
	    this->derivZ = derivZ;

	    d->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);

	    max_xi  = d->gLx;
	    max_eta = d->gLy;
	    max_zta = d->gLz;

	    Nx = d->gNx;
	    Ny = d->gNy;
	    Nz = d->gNz;

	    //Doing this for base y-pencil solvers...
	    c2d->allocY(x);
	    c2d->allocY(y);
	    c2d->allocY(z);

	    //Generate the mesh algebraically...
	    FOR_Z_YPEN{
		FOR_Y_YPEN{
		    FOR_X_YPEN{
			int ip = GETMAJIND_YPEN;
		
			int ii = GETGLOBALXIND_YPEN;
			int jj = GETGLOBALYIND_YPEN;
			int kk = GETGLOBALZIND_YPEN;

			double xi  = d->x[ii];
			double eta = d->y[ii];
			double zta = d->z[ii];
		
			double nXi  = xi/max_xi;
			double nEta = eta/max_eta;
			double nZta = zta/max_zta;

			x[ip] = nXi;
			y[ip] = 2.0*nEta;
			z[ip] = nZta;  

		    }
		}
	    }

	    c2d->allocY(J);
	    c2d->allocY(J11);
	    c2d->allocY(J12);
	    c2d->allocY(J13);
	    c2d->allocY(J21);
	    c2d->allocY(J22);
	    c2d->allocY(J23);
	    c2d->allocY(J31);
	    c2d->allocY(J32);
	    c2d->allocY(J33);

	}


        void solveForJacobians();

};


void AlgebraicSingleBlockMesh::solveForJacobians(){

	double *xE11, *xE21;
	double *xE12, *xE22;
	double *xE13, *xE23;

	c2d->allocY(xE11);
	c2d->allocY(xE12);
	c2d->allocY(xE13);
	c2d->allocY(xE21);
	c2d->allocY(xE22);
	c2d->allocY(xE23);

	//Do the E2 derivatives first...
	derivY->calc1stDerivField(x, xE12);
	derivY->calc1stDerivField(y, xE22);
	
	//Transpose over to E1...
	double *tempX1, *tempX2, *tempX3, *tempX4;
	c2d->allocX(tempX1);
	c2d->allocX(tempX2);
	c2d->allocX(tempX3);
	c2d->allocX(tempX4);

	c2d->transposeY2X_MajorIndex(x, tempX1);
	c2d->transposeY2X_MajorIndex(y, tempX2);

	//Calculate E1 Derivatives..
	derivX->calc1stDerivField(tempX1, tempX3);
	derivX->calc1stDerivField(tempX2, tempX4);

	//Transpose back to E2...
	c2d->transposeX2Y_MajorIndex(tempX3, xE11);
	c2d->transposeX2Y_MajorIndex(tempX4, xE21);

	//Transpose over to E3...
	double *tempZ1, *tempZ2, *tempZ3, *tempZ4;
	c2d->allocZ(tempZ1);
	c2d->allocZ(tempZ2);
	c2d->allocZ(tempZ3);
	c2d->allocZ(tempZ4);

	c2d->transposeY2Z_MajorIndex(x, tempZ1);
	c2d->transposeY2Z_MajorIndex(y, tempZ2);

	//Calculate E3 Derivatives
	derivZ->calc1stDerivField(tempZ1, tempZ3);
	derivZ->calc1stDerivField(tempZ2, tempZ4);

	//Transpose over to E2
	c2d->transposeZ2Y_MajorIndex(tempZ3, xE13);
	c2d->transposeZ2Y_MajorIndex(tempZ4, xE23);


	//Now calculate intermediate components
	double *Va1, *Va2, *Va3;
	double *Vb1, *Vb2, *Vb3;
	double *Vc1, *Vc2, *Vc3;

	c2d->allocY(Va1);
	c2d->allocY(Va2);
	c2d->allocY(Va3);
	c2d->allocY(Vb1);
	c2d->allocY(Vb2);
	c2d->allocY(Vb3);
	c2d->allocY(Vc1);
	c2d->allocY(Vc2);
	c2d->allocY(Vc3);

	FOR_XYZ_YPEN{
	     Va1[ip] = xE11[ip]*y[ip]; 
	     Va2[ip] = xE12[ip]*y[ip]; 
	     Va3[ip] = xE13[ip]*y[ip]; 

	     Vb1[ip] = xE11[ip]*z[ip];
	     Vb2[ip] = xE12[ip]*z[ip];
	     Vb3[ip] = xE13[ip]*z[ip];

	     Vc1[ip] = xE21[ip]*z[ip];
	     Vc2[ip] = xE22[ip]*z[ip];
	     Vc3[ip] = xE23[ip]*z[ip];
	}

	c2d->dealloc(xE11);
	c2d->dealloc(xE12);
	c2d->dealloc(xE13);
	c2d->dealloc(xE21);
	c2d->dealloc(xE22);
	c2d->dealloc(xE23);

	//Only need the off-diagonals of the outer grad tensor of the vectors
	double *dVa12, *dVa13;
	double *dVa21, *dVa23;
	double *dVa31, *dVa32;

	c2d->allocY(dVa12);
	c2d->allocY(dVa13);
	c2d->allocY(dVa21);
	c2d->allocY(dVa23);
	c2d->allocY(dVa31);
	c2d->allocY(dVa32);

	double *dVb12, *dVb13;
	double *dVb21, *dVb23;
	double *dVb31, *dVb32;

	c2d->allocY(dVb12);
	c2d->allocY(dVb13);
	c2d->allocY(dVb21);
	c2d->allocY(dVb23);
	c2d->allocY(dVb31);
	c2d->allocY(dVb32);

	double *dVc12, *dVc13;
	double *dVc21, *dVc23;
	double *dVc31, *dVc32;

	c2d->allocY(dVc12);
	c2d->allocY(dVc13);
	c2d->allocY(dVc21);
	c2d->allocY(dVc23);
	c2d->allocY(dVc31);
	c2d->allocY(dVc32);

	//Start doing the E2 derivatives of all of this stuff
	derivY->calc1stDerivField(Va1, dVa12);
	derivY->calc1stDerivField(Va3, dVa32);
	derivY->calc1stDerivField(Vb1, dVb12);
	derivY->calc1stDerivField(Vb3, dVb32);
	derivY->calc1stDerivField(Vc1, dVc12);
	derivY->calc1stDerivField(Vc3, dVc32);

	//Start doing the E1 derivatives...
	c2d->transposeY2X_MajorIndex(Va2, tempX1);
	c2d->transposeY2X_MajorIndex(Va3, tempX2);
	derivX->calc1stDerivField(tempX1, tempX3);
	derivX->calc1stDerivField(tempX2, tempX4);
	c2d->transposeX2Y_MajorIndex(tempX3, dVa21);
	c2d->transposeX2Y_MajorIndex(tempX4, dVa31);

	c2d->transposeY2X_MajorIndex(Vb2, tempX1);
	c2d->transposeY2X_MajorIndex(Vb3, tempX2);
	derivX->calc1stDerivField(tempX1, tempX3);
	derivX->calc1stDerivField(tempX2, tempX4);
	c2d->transposeX2Y_MajorIndex(tempX3, dVb21);
	c2d->transposeX2Y_MajorIndex(tempX4, dVb31);

	c2d->transposeY2X_MajorIndex(Vc2, tempX1);
	c2d->transposeY2X_MajorIndex(Vc3, tempX2);
	derivX->calc1stDerivField(tempX1, tempX3);
	derivX->calc1stDerivField(tempX2, tempX4);
	c2d->transposeX2Y_MajorIndex(tempX3, dVc21);
	c2d->transposeX2Y_MajorIndex(tempX4, dVc31);

	//Start doing the E3 derivatives...
	c2d->transposeY2Z_MajorIndex(Va1, tempZ1);
	c2d->transposeY2Z_MajorIndex(Va2, tempZ2);
	derivZ->calc1stDerivField(tempZ1, tempZ3);
	derivZ->calc1stDerivField(tempZ2, tempZ4);
	c2d->transposeZ2Y_MajorIndex(tempZ3, dVa13);
	c2d->transposeZ2Y_MajorIndex(tempZ4, dVa23);

	c2d->transposeY2Z_MajorIndex(Vb1, tempZ1);
	c2d->transposeY2Z_MajorIndex(Vb2, tempZ2);
	derivZ->calc1stDerivField(tempZ1, tempZ3);
	derivZ->calc1stDerivField(tempZ2, tempZ4);
	c2d->transposeZ2Y_MajorIndex(tempZ3, dVb13);
	c2d->transposeZ2Y_MajorIndex(tempZ4, dVb23);

	c2d->transposeY2Z_MajorIndex(Vc1, tempZ1);
	c2d->transposeY2Z_MajorIndex(Vc2, tempZ2);
	derivZ->calc1stDerivField(tempZ1, tempZ3);
	derivZ->calc1stDerivField(tempZ2, tempZ4);
	c2d->transposeZ2Y_MajorIndex(tempZ3, dVc13);
	c2d->transposeZ2Y_MajorIndex(tempZ4, dVc23);


	//Start calculating the Jacobian components [(dE/dx)/J]...
	FOR_XYZ_YPEN{
	    J11[ip] = dVc23[ip] - dVc32[ip];
	    J12[ip] = dVb32[ip] - dVb23[ip];
	    J13[ip] = dVa23[ip] - dVa32[ip];

	    J21[ip] = dVc31[ip] - dVc13[ip];
	    J22[ip] = dVb13[ip] - dVb31[ip];
	    J23[ip] = dVa31[ip] - dVa13[ip];

	    J31[ip] = dVc12[ip] - dVc21[ip];
	    J32[ip] = dVb21[ip] - dVb12[ip];
	    J33[ip] = dVa12[ip] - dVa12[ip];
	}

	//Free up all of this space here...
	c2d->deallocXYZ(dVa12);
	c2d->deallocXYZ(dVa13);
	c2d->deallocXYZ(dVa21);
	c2d->deallocXYZ(dVa23);
	c2d->deallocXYZ(dVa31);
	c2d->deallocXYZ(dVa32);

	c2d->deallocXYZ(dVb12);
	c2d->deallocXYZ(dVb13);
	c2d->deallocXYZ(dVb21);
	c2d->deallocXYZ(dVb23);
	c2d->deallocXYZ(dVb31);
	c2d->deallocXYZ(dVb32);

	c2d->deallocXYZ(dVc12);
	c2d->deallocXYZ(dVc13);
	c2d->deallocXYZ(dVc21);
	c2d->deallocXYZ(dVc23);
	c2d->deallocXYZ(dVc31);
	c2d->deallocXYZ(dVc32);

	c2d->deallocXYZ(tempX3);
	c2d->deallocXYZ(tempX4);
	c2d->deallocXYZ(tempZ3);
	c2d->deallocXYZ(tempZ4);


	//Now start computing the Jacobian determinant...
	double *preJdet1, *preJdet2, *preJdet3;
	double *Jdet1, *Jdet2, *Jdet3;
	c2d->allocY(preJdet1);
	c2d->allocY(preJdet2);
	c2d->allocY(preJdet3);
	c2d->allocY(Jdet1);
	c2d->allocY(Jdet2);
	c2d->allocY(Jdet3);

	FOR_XYZ_YPEN{
	    preJdet1[ip] = x[ip]*J11[ip] + y[ip]*J12[ip] + z[ip]*J13[ip];
	    preJdet2[ip] = x[ip]*J21[ip] + y[ip]*J22[ip] + z[ip]*J23[ip];
	    preJdet3[ip] = x[ip]*J31[ip] + y[ip]*J32[ip] + z[ip]*J33[ip];
	}

	//Compute the E2 component...
	derivY->calc1stDerivField(preJdet2, Jdet2);

	//Compute the E1 component....
	c2d->transposeY2X_MajorIndex(preJdet1, tempX1);
	derivX->calc1stDerivField(tempX1, tempX2);
	c2d->transposeX2Y_MajorIndex(tempX2, Jdet1);

 	//Compute the E3 component....
	c2d->transposeY2Z_MajorIndex(preJdet3, tempZ1);
	derivZ->calc1stDerivField(tempZ1, tempZ2);
	c2d->transposeZ2Y_MajorIndex(tempZ2, Jdet3);

	//Compute the Jacobian derivative...
	FOR_XYZ_YPEN{
	    J[ip] = (1.0/3.0)*(Jdet1[ip] + Jdet2[ip] + Jdet3[ip]);
	}

	//Free up all of the spaces we've been using...
	c2d->deallocXYZ(tempX1);
	c2d->deallocXYZ(tempX2);
	c2d->deallocXYZ(tempZ1);
	c2d->deallocXYZ(tempZ2);

	c2d->deallocXYZ(preJdet1);
	c2d->deallocXYZ(preJdet2);
	c2d->deallocXYZ(preJdet3);

	c2d->deallocXYZ(Jdet1);
	c2d->deallocXYZ(Jdet2);
	c2d->deallocXYZ(Jdet3);

}

#endif
