#ifndef _CALGEBRAICSINGLEBLOCKMESHH_
#define _CALGEBRAICSINGLEBLOCKMESHH_

#include "Macros.hpp"
#include "Utils.hpp"
#include "AbstractSingleBlockMesh.hpp"

class AlgebraicSingleBlockMesh:public AbstractSingleBlockMesh{

    public:

        C2Decomp *c2d;
        AbstractCSolver *cs;
        Domain *d;

	AlgebraicSingleBlockMesh(C2Decomp *c2d, AbstractCSolver *cs, Domain *dom, int mpiRank){


	    this->mpiRank = mpiRank;
	    this->cs = cs;
	    this->c2d = c2d;
	    this->d = dom;
	    this->derivX = cs->derivX;
	    this->derivY = cs->derivY;
	    this->derivZ = cs->derivZ;

	    d->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);

	    max_xi  = d->gLx;
	    max_eta = d->gLy;
	    max_zta = d->gLz;

	    Nx = d->gNx;
	    Ny = d->gNy;
	    Nz = d->gNz;


	    IF_RANK0 cout << endl << " > Running Algebraic single block mesh initialization for curvlinear grid solver! " << endl;

	    //Doing this for base y-pencil solvers...
	    c2d->allocY(x);
	    c2d->allocY(y);
	    c2d->allocY(z);


	    //Get the global mins/maxes
	    double x_min_local[3] = {1e10, 1e10, 1e10};
	    double x_max_local[3] = {-1e10, -1e10, -1e10};

	    //Generate the mesh algebraically...
	    FOR_Z_YPEN{
		FOR_Y_YPEN{
		    FOR_X_YPEN{
			int ip = GETMAJIND_YPEN;
		
			int ii = GETGLOBALXIND_YPEN;
			int jj = GETGLOBALYIND_YPEN;
			int kk = GETGLOBALZIND_YPEN;

			double xi  = d->x[ii];
			double eta = d->y[jj];
			double zta = d->z[kk];
		
			//double nXi  = xi/max_xi;
			//double nEta = eta/max_eta;
			//double nZta = zta/max_zta;

			x[ip] = xi+1.5*eta;
			y[ip] = eta;
			z[ip] = zta;  

			//Since we're already in this loop, calculate the local max and mins
			x_min_local[0] = fmin(x_min_local[0], x[ip]);
			x_min_local[1] = fmin(x_min_local[1], y[ip]);
			x_min_local[2] = fmin(x_min_local[2], z[ip]);
	
			x_max_local[0] = fmax(x_max_local[0], x[ip]);
			x_max_local[1] = fmax(x_max_local[1], y[ip]);
			x_max_local[2] = fmax(x_max_local[2], z[ip]);

		    }
		}
	    }

	    //Reduce to get the global bounding box
	    MPI_Allreduce(x_min_local, x_min, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	    MPI_Allreduce(x_max_local, x_max, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);


	    if(cs->bc->bcXType == BC::PERIODIC_SOLVE){
		periodicX = true;
		periodicXTranslation[0] = 1.0;
		periodicXTranslation[1] = 0.0;
		periodicXTranslation[2] = 0.0;
		IF_RANK0 cout << "  Periodic x-face translation = {" << periodicXTranslation[0] << ", " << periodicXTranslation[1] << ", " << periodicXTranslation[2] << "}" << endl;;

	    }else{
		periodicX = false;
	    }

	    if(cs->bc->bcYType == BC::PERIODIC_SOLVE){
		periodicY = true;
		periodicYTranslation[0] = 1.5;
		periodicYTranslation[1] = 1.0;
		periodicYTranslation[2] = 0.0;
		IF_RANK0 cout << "  Periodic y-face translation = {" << periodicYTranslation[0] << ", " << periodicYTranslation[1] << ", " << periodicYTranslation[2] << "}" << endl;

	    }else{
		periodicY = false; 
	    }

	    if(cs->bc->bcZType == BC::PERIODIC_SOLVE){
		periodicZ = true;
		periodicZTranslation[0] = 0.0;
		periodicZTranslation[1] = 0.0;
		periodicZTranslation[2] = 1.0;

		IF_RANK0 cout << "  Periodic z-face translation = {" << periodicZTranslation[0] << ", " << periodicZTranslation[1] << ", " << periodicZTranslation[2] << "}" << endl;;

	    }else{
		periodicZ = false;
	    }

	    getRange(x, "X", pySize[0], pySize[1], pySize[2], mpiRank);
	    getRange(y, "Y", pySize[0], pySize[1], pySize[2], mpiRank);
	    getRange(z, "Z", pySize[0], pySize[1], pySize[2], mpiRank);
	    IF_RANK0 cout << "  Range of Xi: "  << d->x[0] <<  ":" << d->x[Nx-1] << endl;
	    IF_RANK0 cout << "  Range of Eta: " << d->y[0] <<  ":" << d->y[Ny-1] << endl;
	    IF_RANK0 cout << "  Range of Zta: " << d->z[0] <<  ":" << d->z[Nz-1] << endl;


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


	    //Initialize the ADT object for interpolation and locating points in the grid...
	    IF_RANK0 cout << " > Initializing ADT..." << endl;

            //Get our coordinates in neighbors across partitions using halo updates
            double *x_halo = NULL;
            double *y_halo = NULL;
            double *z_halo = NULL;

	    generateCoordinateHaloArrays(x_halo, y_halo, z_halo);
	    
	    int Nlocal = pySize[0]*pySize[1]*pySize[2];
            double (*boundBoxMin)[3] = new double[Nlocal][3];
            double (*boundBoxMax)[3] = new double[Nlocal][3];

            //Cycle through the halo arrays of coordinates
            for(int kp = 0; kp < pySize[2]; kp++){
                for(int jp = 0; jp < pySize[1]; jp++){
                    for(int ip = 0; ip < pySize[0]; ip++){

                        //This is the non-halo array index
                        int ii = kp*pySize[1]*pySize[0] + jp*pySize[0] + ip;

                        //This is the non-halo array index, y-major
                        int ii_major = ip*pySize[2]*pySize[1] + kp*pySize[1] + jp;


			double box_p[8][3];
			getOrderedBlockCoordinates(ip, jp, kp, x_halo, y_halo, z_halo, box_p);

			double x_max = -1.0e100; 
			double y_max = -1.0e100; 
			double z_max = -1.0e100; 
			double x_min =  1.0e100;
			double y_min =  1.0e100;
			double z_min =  1.0e100;


		        bool xEndFlag = false;
			if(pyStart[0] + ip == Nx-1){
			    xEndFlag = true;
			}

			bool yEndFlag = false;
			if(jp == Ny-1){
			    yEndFlag = true;
			}

			bool zEndFlag = false;
			if(pyStart[2] + kp == Nz-1){
			    zEndFlag = true;
			}      

			bool noBBFlag = false;	
                        if((xEndFlag && !periodicX) ||
			   (yEndFlag && !periodicY) ||
			   (zEndFlag && !periodicZ)){
                            noBBFlag = true;
                        }

			if(!noBBFlag){

			    for(int iip = 0; iip < 8; iip++){
    			       x_max = fmax(x_max, box_p[iip][0]);
    			       x_min = fmin(x_min, box_p[iip][0]);

    			       y_max = fmax(y_max, box_p[iip][1]);
    			       y_min = fmin(y_min, box_p[iip][1]);

    			       z_max = fmax(z_max, box_p[iip][2]);
    			       z_min = fmin(z_min, box_p[iip][2]);
			    }

		  	}else{
	 		    x_max = x[ii];
			    x_min = x[ii];
			    y_max = y[ii];
			    y_min = y[ii];
			    z_min = z[ii];
			    z_min = z[ii];
			}

			//We'll usually be accessing this in the major indexing fashion
             	        boundBoxMin[ii_major][0] = x_min;
             	        boundBoxMin[ii_major][1] = y_min;
             	        boundBoxMin[ii_major][2] = z_min;
			
;
		        boundBoxMax[ii_major][0] = x_max; 
		        boundBoxMax[ii_major][1] = y_max; 
		        boundBoxMax[ii_major][2] = z_max; 
                    }
                }
            }

	    FOR_XYZ_YPEN{
		FOR_I3{
		    double delta = 1.0E-6*(boundBoxMax[ip][i] - boundBoxMin[ip][i]);
		    boundBoxMax[ip][i] += delta;
		    boundBoxMin[ip][i] -= delta;
		}
	    }

	    IF_RANK0 cout << " > Done getting bounding boxes for the CV's, initializing ADT... " << endl;

	    adt = new Adt<double>(Nlocal, boundBoxMin, boundBoxMax);

	    IF_RANK0 cout << " > Done!" << endl;


	    delete[] boundBoxMin;
	    delete[] boundBoxMax;
	    delete[] x_halo;
	    delete[] y_halo;
	    delete[] z_halo;

	}

 	void solveForJacobians();

	void getOrderedBlockCoordinates(int ip, int jp, int kp, double *x_halo, double *y_halo, double *z_halo, double box_p[8][3]);
	void getOrderedBlockXiCoordinates(int ip, int jp, int kp, double box_pxi[8][3]);

	int findCVForPoint(double p[3], double *x_halo, double *y_halo, double *z_halo);

	void generateCoordinateHaloArrays(double *&x_halo, double *&y_halo, double *&z_halo);

};


void AlgebraicSingleBlockMesh::getOrderedBlockXiCoordinates(int ip, int jp, int kp, double box_pxi[8][3]){

        double dxi, deta, dzta;

        //This is correct if the domain bounds have been properly set to their max xi, eta, and zeta values in the
        //MPI_Compact.cpp file
        dxi  = d->dx;
        deta = d->dy;
        dzta = d->dz;

        box_pxi[0][0] = dxi *(double)ip;
        box_pxi[1][0] = dxi *(double)ip;
        box_pxi[2][0] = dxi *(double)ip;
        box_pxi[3][0] = dxi *(double)ip;
        box_pxi[4][0] = dxi *(double)(ip+1);
        box_pxi[5][0] = dxi *(double)(ip+1);
        box_pxi[6][0] = dxi *(double)(ip+1);
        box_pxi[7][0] = dxi *(double)(ip+1);

        box_pxi[0][1] = deta*(double)jp;
        box_pxi[1][1] = deta*(double)jp;
        box_pxi[2][1] = deta*(double)(jp+1);
        box_pxi[3][1] = deta*(double)(jp+1);
        box_pxi[4][1] = deta*(double)jp;
        box_pxi[5][1] = deta*(double)jp;
        box_pxi[6][1] = deta*(double)(jp+1);
        box_pxi[7][1] = deta*(double)(jp+1);

        box_pxi[0][2] = dzta*(double)kp;
        box_pxi[1][2] = dzta*(double)(kp+1);
        box_pxi[2][2] = dzta*(double)kp;
        box_pxi[3][2] = dzta*(double)(kp+1);
        box_pxi[4][2] = dzta*(double)kp;
        box_pxi[5][2] = dzta*(double)(kp+1);
        box_pxi[6][2] = dzta*(double)kp;
        box_pxi[7][2] = dzta*(double)(kp+1);

};

void AlgebraicSingleBlockMesh::generateCoordinateHaloArrays(double *&x_halo, double *&y_halo, double *&z_halo){

	//Really should implement halo transfers for major indexed arrays
        double *x_temp1, *y_temp1, *z_temp1;
        c2d->allocX(x_temp1);
        c2d->allocX(y_temp1);
        c2d->allocX(z_temp1);

        //So move to x-pencil then to non y-major y_pencil...
        c2d->transposeY2X_MajorIndex(x, x_temp1);
        c2d->transposeY2X_MajorIndex(y, y_temp1);
        c2d->transposeY2X_MajorIndex(z, z_temp1);

        //Then back over to y-pencil in x-major array...
        double *x_temp2, *y_temp2, *z_temp2;
        c2d->allocY(x_temp2); c2d->allocY(y_temp2); c2d->allocY(z_temp2);

        c2d->transposeX2Y(x_temp1, x_temp2);
        c2d->transposeX2Y(y_temp1, y_temp2);
        c2d->transposeX2Y(z_temp1, z_temp2);

        delete[] x_temp1;
        delete[] y_temp1;
        delete[] z_temp1;

        c2d->updateHalo(x_temp2, x_halo, 1, 1);
        c2d->updateHalo(y_temp2, y_halo, 1, 1);
        c2d->updateHalo(z_temp2, z_halo, 1, 1);

        delete[] x_temp2;
        delete[] y_temp2;
        delete[] z_temp2;
};

void AlgebraicSingleBlockMesh::getOrderedBlockCoordinates(int ip, int jp, int kp, double *x_halo, double *y_halo, double *z_halo, double box_p[8][3]){

	int iih_0_0_0;
	int iih_0_0_1;
	int iih_0_1_0;
	int iih_0_1_1;
	int iih_1_0_0;
	int iih_1_0_1;
	int iih_1_1_0;
	int iih_1_1_1;

	//What if we're trying to access *p+1 and we're not periodic in that direction? what happens now in x and z? What should 
	//we return for y? just the origin coordinate?
	//No its fine because C2Decomp will allocate the array with padding even if its not periodic and the *p+1 point should never
	//be returned since we've zero'd the bounding box from above


        //This is the halo array index for the same point
        iih_0_0_0 = (kp+1)*pySize[1]*(pySize[0]+2) + jp*(pySize[0]+2) + ip + 1;

        //Halo array index for i, j, k+1
        iih_0_0_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (jp)*(pySize[0]+2) + ip + 1;

        //Halo array index for i, j+1, k
	if(periodicY && jp == (Ny-1)){
            iih_0_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (0)*(pySize[0]+2) + ip + 1;
	}else{
            iih_0_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 1;
	}

        //Halo array index for i, j+1, k+1
	if(periodicY && jp == (Ny-1)){
            iih_0_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (0)*(pySize[0]+2) + ip + 1;
	}else{
            iih_0_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 1;
	}

        //Halo array index for i+1, j, k
        iih_1_0_0 = (kp+1)*pySize[1]*(pySize[0]+2) + jp*(pySize[0]+2) + ip + 2;

        //Halo array index for i+1, j, k+1
        iih_1_0_1 = (kp+2)*pySize[1]*(pySize[0]+2) + jp*(pySize[0]+2) + ip + 2;

        //Halo array index for i+1, j+1, k
	if(periodicY && jp == (Ny-1)){
            iih_1_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (0)*(pySize[0]+2) + ip + 2;
	}else{
            iih_1_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 2;
	}

        //Halo array index for i+1, j+1, k+1
	if(periodicY && jp == (Ny-1)){
	    iih_1_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (0)*(pySize[0]+2) + ip + 2;
	}else{
	    iih_1_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 2;
	}
 
        bool xEndFlag = false;
        if(pyStart[0] + ip == Nx-1){
            xEndFlag = true;
        }

        bool yEndFlag = false;
        if(jp == Ny-1){
            yEndFlag = true;
        }

        bool zEndFlag = false;
        if(pyStart[2] + kp == Nz-1){
            zEndFlag = true;
        }      

	/////////////////////////////////
	//Definite local point first...//
	//Do 0 0 0 first/////////////////
	//////////////////
	box_p[0][0] = x_halo[iih_0_0_0];
	box_p[0][1] = y_halo[iih_0_0_0];
	box_p[0][2] = z_halo[iih_0_0_0];


	/////////////////
	//Next do 0 0 1// 
	/////////////////

	if(zEndFlag && periodicZ){
	    box_p[1][0] = x_halo[iih_0_0_1] + periodicZTranslation[0];
	    box_p[1][1] = y_halo[iih_0_0_1] + periodicZTranslation[1];
	    box_p[1][2] = z_halo[iih_0_0_1] + periodicZTranslation[2];
	}else{ //If here, we're an interior point
	    box_p[1][0] = x_halo[iih_0_0_1];
	    box_p[1][1] = y_halo[iih_0_0_1];
	    box_p[1][2] = z_halo[iih_0_0_1];
	}



	/////////////////
	//Next do 0 1 0//
	/////////////////

	if(yEndFlag && periodicY){
	    box_p[2][0] = x_halo[iih_0_1_0] + periodicYTranslation[0];
	    box_p[2][1] = y_halo[iih_0_1_0] + periodicYTranslation[1];
	    box_p[2][2] = z_halo[iih_0_1_0] + periodicYTranslation[2];
	}else{// If, here we're an interior point
	    box_p[2][0] = x_halo[iih_0_1_0];
	    box_p[2][1] = y_halo[iih_0_1_0];
	    box_p[2][2] = z_halo[iih_0_1_0];
	}

	/////////////////
	//Next do 0 1 1//
	/////////////////

	if(zEndFlag && periodicZ){
	    if(yEndFlag && periodicY){
	        box_p[3][0] = x_halo[iih_0_1_1] + periodicZTranslation[0] + periodicYTranslation[0];
	        box_p[3][1] = y_halo[iih_0_1_1] + periodicZTranslation[1] + periodicYTranslation[1];
	        box_p[3][2] = z_halo[iih_0_1_1] + periodicZTranslation[2] + periodicYTranslation[2];
	    }else{
		box_p[3][0] = x_halo[iih_0_1_1] + periodicZTranslation[0];
		box_p[3][1] = y_halo[iih_0_1_1] + periodicZTranslation[1];
		box_p[3][2] = z_halo[iih_0_1_1] + periodicZTranslation[2];
	    }
	}else{ // in interior domain in z-direction
	    if(yEndFlag && periodicY){
		box_p[3][0] = x_halo[iih_0_1_1] + periodicYTranslation[0];
		box_p[3][1] = y_halo[iih_0_1_1] + periodicYTranslation[1];
		box_p[3][2] = z_halo[iih_0_1_1] + periodicYTranslation[2];
	    }else{//If we're here, we're interior 
		box_p[3][0] = x_halo[iih_0_1_1];
		box_p[3][1] = y_halo[iih_0_1_1];
		box_p[3][2] = z_halo[iih_0_1_1];
	    }
	}


	/////////////////
	//Next do 1 0 0//
	/////////////////

	if(xEndFlag && periodicX){
	    box_p[4][0] = x_halo[iih_1_0_0] + periodicXTranslation[0];
	    box_p[4][1] = y_halo[iih_1_0_0] + periodicXTranslation[1];
	    box_p[4][2] = z_halo[iih_1_0_0] + periodicXTranslation[2];
	}else{
	    box_p[4][0] = x_halo[iih_1_0_0];
	    box_p[4][1] = y_halo[iih_1_0_0];
	    box_p[4][2] = z_halo[iih_1_0_0];
	}

	/////////////////
	//Next do 1 0 1//
	/////////////////

	if(xEndFlag && periodicX){
	    if(zEndFlag && periodicZ){
	        box_p[5][0] = x_halo[iih_1_0_1] + periodicXTranslation[0] + periodicZTranslation[0];
	        box_p[5][1] = y_halo[iih_1_0_1] + periodicXTranslation[1] + periodicZTranslation[1];
	        box_p[5][2] = z_halo[iih_1_0_1] + periodicXTranslation[2] + periodicZTranslation[2];
	    }else{
		box_p[5][0] = x_halo[iih_1_0_1] + periodicXTranslation[0];
		box_p[5][1] = y_halo[iih_1_0_1] + periodicXTranslation[1];
		box_p[5][2] = z_halo[iih_1_0_1] + periodicXTranslation[2];
	    }
	}else{
	    if(zEndFlag && periodicZ){
		box_p[5][0] = x_halo[iih_1_0_1] + periodicZTranslation[0];
		box_p[5][1] = y_halo[iih_1_0_1] + periodicZTranslation[1];
		box_p[5][2] = z_halo[iih_1_0_1] + periodicZTranslation[2];
	    }else{
		box_p[5][0] = x_halo[iih_1_0_1];
		box_p[5][1] = y_halo[iih_1_0_1];
		box_p[5][2] = z_halo[iih_1_0_1];
	    }
	}

        /////////////////
	//Next do 1 1 0// 
	/////////////////   

	if(xEndFlag && periodicX){
	    if(yEndFlag && periodicY){
		box_p[6][0] = x_halo[iih_1_1_0] + periodicXTranslation[0] + periodicYTranslation[0];
		box_p[6][1] = y_halo[iih_1_1_0] + periodicXTranslation[1] + periodicYTranslation[1];
		box_p[6][2] = z_halo[iih_1_1_0] + periodicXTranslation[2] + periodicYTranslation[2];
	    }else{
		box_p[6][0] = x_halo[iih_1_1_0] + periodicXTranslation[0];
		box_p[6][1] = y_halo[iih_1_1_0] + periodicXTranslation[1];
		box_p[6][2] = z_halo[iih_1_1_0] + periodicXTranslation[2];
	    }
	}else{
	    if(yEndFlag && periodicY){
		box_p[6][0] = x_halo[iih_1_1_0] + periodicYTranslation[0];
		box_p[6][1] = y_halo[iih_1_1_0] + periodicYTranslation[1];
		box_p[6][2] = z_halo[iih_1_1_0] + periodicYTranslation[2];
	    }else{
		box_p[6][0] = x_halo[iih_1_1_0];
		box_p[6][1] = y_halo[iih_1_1_0];
		box_p[6][2] = z_halo[iih_1_1_0];
	    }
	}


        ////////////////////
        //Finally do 1 1 1//
        ////////////////////

        if(xEndFlag && periodicX){
	    if(yEndFlag && periodicY){
	        if(zEndFlag && periodicZ){
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicXTranslation[0] + periodicYTranslation[0] + periodicZTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicXTranslation[1] + periodicYTranslation[1] + periodicZTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicXTranslation[2] + periodicYTranslation[2] + periodicZTranslation[2];
	        }else{
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicXTranslation[0] + periodicYTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicXTranslation[1] + periodicYTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicXTranslation[2] + periodicYTranslation[2];
	        }
	    }else{
	        if(zEndFlag && periodicZ){
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicXTranslation[0] + periodicZTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicXTranslation[1] + periodicZTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicXTranslation[2] + periodicZTranslation[2];
	        }else{
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicXTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicXTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicXTranslation[2];
	        }
	    }
        }else{
	    if(yEndFlag && periodicY){
	        if(zEndFlag && periodicZ){
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicYTranslation[0] + periodicZTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicYTranslation[1] + periodicZTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicYTranslation[2] + periodicZTranslation[2];
	        }else{
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicYTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicYTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicYTranslation[2];
	        }
	    }else{
	        if(zEndFlag && periodicZ){
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicZTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicZTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicZTranslation[2];
	        }else{
		    box_p[7][0] = x_halo[iih_1_1_1];
		    box_p[7][1] = y_halo[iih_1_1_1];
		    box_p[7][2] = z_halo[iih_1_1_1];
	        }
	    }
        }


};

int AlgebraicSingleBlockMesh::findCVForPoint(double p[3], double *x_halo, double *y_halo, double *z_halo){

    int cvListSize, cvList[8192];

    adt->buildListForPoint(cvListSize, cvList, p);

    //This function will return -1 for a point not found and the local control volume index in y-major
    //indexing if it is found
    int base_index = -1;
    for(int ii = 0; ii < cvListSize; ii++){

	base_index = cvList[ii];

	//Back out the ip, jp, kp coordinates based off single major index int
	int jp =  base_index%pySize[1];
	int kp = (base_index/pySize[1])%pySize[2];
	int ip =  base_index/(pySize[2]*pySize[1]);

	double box_p[8][3];

	getOrderedBlockCoordinates(ip, jp, kp, x_halo, y_halo, z_halo, box_p);	

	if(isPointInHexa(p, box_p)){
	    base_index = cvList[ii];
	    break;
	}else{
	    base_index = -1;
	}
    }

    return base_index;    

}


void AlgebraicSingleBlockMesh::solveForJacobians(){


	IF_RANK0 cout << " > Solving for Jacobian Matrix... ";

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
	if(periodicY){
	    double *Nm2x, *Nm1x, *Np1x, *Np2x;
	    Nm2x = new double[pySize[0]*pySize[2]];
	    Nm1x = new double[pySize[0]*pySize[2]];
	    Np1x = new double[pySize[0]*pySize[2]];
	    Np2x = new double[pySize[0]*pySize[2]];

	    double *Nm2y, *Nm1y, *Np1y, *Np2y;
	    Nm2y = new double[pySize[0]*pySize[2]];
	    Nm1y = new double[pySize[0]*pySize[2]];
	    Np1y = new double[pySize[0]*pySize[2]];
	    Np2y = new double[pySize[0]*pySize[2]];

	    FOR_X_YPEN{
		FOR_Z_YPEN{
		    int ii = i*pySize[2] + k;

		    int iim2 = i*pySize[2]*pySize[1] + k*pySize[1] + pySize[1]-2;		
		    int iim1 = i*pySize[2]*pySize[1] + k*pySize[1] + pySize[1]-1;		
		    int iip1 = i*pySize[2]*pySize[1] + k*pySize[1] + 0;		
		    int iip2 = i*pySize[2]*pySize[1] + k*pySize[1] + 1;		

		    Nm2x[ii] = x[iim2]-periodicYTranslation[0];
		    Nm1x[ii] = x[iim1]-periodicYTranslation[0];
		    Np1x[ii] = x[iip1]+periodicYTranslation[0];
		    Np2x[ii] = x[iip2]+periodicYTranslation[0];
	
		    Nm2y[ii] = y[iim2]-periodicYTranslation[1];
		    Nm1y[ii] = y[iim1]-periodicYTranslation[1];
		    Np1y[ii] = y[iip1]+periodicYTranslation[1];
		    Np2y[ii] = y[iip2]+periodicYTranslation[1];
	
		}
	    }

	    derivY->calc1stDerivField_TPB(x, xE12, Nm2x, Nm1x, Np1x, Np2x);
	    derivY->calc1stDerivField_TPB(y, xE22, Nm2y, Nm1y, Np1y, Np2y);

	    delete[] Nm2x;
	    delete[] Nm1x;
	    delete[] Np1x;
	    delete[] Np2x;

	    delete[] Nm2y;
	    delete[] Nm1y;
	    delete[] Np1y;
	    delete[] Np2y;

	}else{
	    derivY->calc1stDerivField(x, xE12);
	    derivY->calc1stDerivField(y, xE22);
	}
	

	//Transpose over to E1...
	double *tempX1, *tempX2, *tempX3, *tempX4;
	c2d->allocX(tempX1);
	c2d->allocX(tempX2);
	c2d->allocX(tempX3);
	c2d->allocX(tempX4);

	c2d->transposeY2X_MajorIndex(x, tempX1);
	c2d->transposeY2X_MajorIndex(y, tempX2);

	//Calculate E1 Derivatives..
	if(periodicX){
	    double *Nm2x, *Nm1x, *Np1x, *Np2x;
            Nm2x = new double[pxSize[1]*pxSize[2]];
            Nm1x = new double[pxSize[1]*pxSize[2]];
            Np1x = new double[pxSize[1]*pxSize[2]];
            Np2x = new double[pxSize[1]*pxSize[2]];

	    double *Nm2y, *Nm1y, *Np1y, *Np2y;
            Nm2y = new double[pxSize[1]*pxSize[2]];
            Nm1y = new double[pxSize[1]*pxSize[2]];
            Np1y = new double[pxSize[1]*pxSize[2]];
            Np2y = new double[pxSize[1]*pxSize[2]];

            FOR_Z_XPEN{
                FOR_Y_XPEN{
                    int ii = k*pxSize[1] + j;

                    int iim2 = k*pxSize[0]*pxSize[1] + j*pxSize[0] + pxSize[0]-2;
                    int iim1 = k*pxSize[0]*pxSize[1] + j*pxSize[0] + pxSize[0]-1;
                    int iip1 = k*pxSize[0]*pxSize[1] + j*pxSize[0] + 0;
                    int iip2 = k*pxSize[0]*pxSize[1] + j*pxSize[0] + 1;

                    Nm2x[ii] = tempX1[iim2]-periodicXTranslation[0];
                    Nm1x[ii] = tempX1[iim1]-periodicXTranslation[0];
                    Np1x[ii] = tempX1[iip1]+periodicXTranslation[0];
                    Np2x[ii] = tempX1[iip2]+periodicXTranslation[0];

                    Nm2y[ii] = tempX2[iim2]-periodicXTranslation[1];
                    Nm1y[ii] = tempX2[iim1]-periodicXTranslation[1];
                    Np1y[ii] = tempX2[iip1]+periodicXTranslation[1];
                    Np2y[ii] = tempX2[iip2]+periodicXTranslation[1];

                }
            }

            derivX->calc1stDerivField_TPB(tempX1, tempX3, Nm2x, Nm1x, Np1x, Np2x);
            derivX->calc1stDerivField_TPB(tempX2, tempX4, Nm2y, Nm1y, Np1y, Np2y);

            delete[] Nm2x;
            delete[] Nm1x;
            delete[] Np1x;
            delete[] Np2x;

            delete[] Nm2y;
            delete[] Nm1y;
            delete[] Np1y;
            delete[] Np2y;

	}else{
	    derivX->calc1stDerivField(tempX1, tempX3);
	    derivX->calc1stDerivField(tempX2, tempX4);
	}


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
	if(periodicZ){
	    double *Nm2x, *Nm1x, *Np1x, *Np2x;
            Nm2x = new double[pzSize[1]*pzSize[0]];
            Nm1x = new double[pzSize[1]*pzSize[0]];
            Np1x = new double[pzSize[1]*pzSize[0]];
            Np2x = new double[pzSize[1]*pzSize[0]];

	    double *Nm2y, *Nm1y, *Np1y, *Np2y;
            Nm2y = new double[pzSize[1]*pzSize[0]];
            Nm1y = new double[pzSize[1]*pzSize[0]];
            Np1y = new double[pzSize[1]*pzSize[0]];
            Np2y = new double[pzSize[1]*pzSize[0]];

            FOR_Y_ZPEN{
                FOR_X_ZPEN{
                    int ii = j*pzSize[0] + i;

                    int iim2 = j*pzSize[0]*pzSize[2] + i*pzSize[2] + pzSize[2]-2;
                    int iim1 = j*pzSize[0]*pzSize[2] + i*pzSize[2] + pzSize[2]-1;
                    int iip1 = j*pzSize[0]*pzSize[2] + i*pzSize[2] + 0;
                    int iip2 = j*pzSize[0]*pzSize[2] + i*pzSize[2] + 1;

                    Nm2x[ii] = tempZ1[iim2]-periodicZTranslation[0];
                    Nm1x[ii] = tempZ1[iim1]-periodicZTranslation[0];
                    Np1x[ii] = tempZ1[iip1]+periodicZTranslation[0];
                    Np2x[ii] = tempZ1[iip2]+periodicZTranslation[0];

                    Nm2y[ii] = tempZ2[iim2]-periodicZTranslation[1];
                    Nm1y[ii] = tempZ2[iim1]-periodicZTranslation[1];
                    Np1y[ii] = tempZ2[iip1]+periodicZTranslation[1];
                    Np2y[ii] = tempZ2[iip2]+periodicZTranslation[1];

                }
            }

            derivZ->calc1stDerivField_TPB(tempZ1, tempZ3, Nm2x, Nm1x, Np1x, Np2x);
            derivZ->calc1stDerivField_TPB(tempZ2, tempZ4, Nm2y, Nm1y, Np1y, Np2y);

            delete[] Nm2x;
            delete[] Nm1x;
            delete[] Np1x;
            delete[] Np2x;

            delete[] Nm2y;
            delete[] Nm1y;
            delete[] Np1y;
            delete[] Np2y;


	}else{
	    derivZ->calc1stDerivField(tempZ1, tempZ3);
	    derivZ->calc1stDerivField(tempZ2, tempZ4);
	}

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
	if(periodicY){
	    double *Nm2a1, *Nm1a1, *Np1a1, *Np2a1;
	    Nm2a1 = new double[pySize[0]*pySize[2]];
	    Nm1a1 = new double[pySize[0]*pySize[2]];
	    Np1a1 = new double[pySize[0]*pySize[2]];
	    Np2a1 = new double[pySize[0]*pySize[2]];

	    double *Nm2a3, *Nm1a3, *Np1a3, *Np2a3;
	    Nm2a3 = new double[pySize[0]*pySize[2]];
	    Nm1a3 = new double[pySize[0]*pySize[2]];
	    Np1a3 = new double[pySize[0]*pySize[2]];
	    Np2a3 = new double[pySize[0]*pySize[2]];

	    double *Nm2b1, *Nm1b1, *Np1b1, *Np2b1;
	    Nm2b1 = new double[pySize[0]*pySize[2]];
	    Nm1b1 = new double[pySize[0]*pySize[2]];
	    Np1b1 = new double[pySize[0]*pySize[2]];
	    Np2b1 = new double[pySize[0]*pySize[2]];

	    double *Nm2b3, *Nm1b3, *Np1b3, *Np2b3;
	    Nm2b3 = new double[pySize[0]*pySize[2]];
	    Nm1b3 = new double[pySize[0]*pySize[2]];
	    Np1b3 = new double[pySize[0]*pySize[2]];
	    Np2b3 = new double[pySize[0]*pySize[2]];

	    double *Nm2c1, *Nm1c1, *Np1c1, *Np2c1;
	    Nm2c1 = new double[pySize[0]*pySize[2]];
	    Nm1c1 = new double[pySize[0]*pySize[2]];
	    Np1c1 = new double[pySize[0]*pySize[2]];
	    Np2c1 = new double[pySize[0]*pySize[2]];

	    double *Nm2c3, *Nm1c3, *Np1c3, *Np2c3;
	    Nm2c3 = new double[pySize[0]*pySize[2]];
	    Nm1c3 = new double[pySize[0]*pySize[2]];
	    Np1c3 = new double[pySize[0]*pySize[2]];
	    Np2c3 = new double[pySize[0]*pySize[2]];


	    FOR_X_YPEN{
		FOR_Z_YPEN{
		    int ii = i*pySize[2] + k;

		    int iim2 = i*pySize[2]*pySize[1] + k*pySize[1] + pySize[1]-2;		
		    int iim1 = i*pySize[2]*pySize[1] + k*pySize[1] + pySize[1]-1;		
		    int iip1 = i*pySize[2]*pySize[1] + k*pySize[1] + 0;		
		    int iip2 = i*pySize[2]*pySize[1] + k*pySize[1] + 1;		

		    Nm2a1[ii] = (y[iim2]-periodicYTranslation[1])*xE11[iim2];
		    Nm1a1[ii] = (y[iim1]-periodicYTranslation[1])*xE11[iim1];
		    Np1a1[ii] = (y[iip1]+periodicYTranslation[1])*xE11[iip1];
		    Np2a1[ii] = (y[iip2]+periodicYTranslation[1])*xE11[iip2];
	
		    Nm2a3[ii] = (y[iim2]-periodicYTranslation[1])*xE13[iim2];
		    Nm1a3[ii] = (y[iim1]-periodicYTranslation[1])*xE13[iim1];
		    Np1a3[ii] = (y[iip1]+periodicYTranslation[1])*xE13[iip1];
		    Np2a3[ii] = (y[iip2]+periodicYTranslation[1])*xE13[iip2];

		    Nm2b1[ii] = (z[iim2]-periodicYTranslation[2])*xE11[iim2];
		    Nm1b1[ii] = (z[iim1]-periodicYTranslation[2])*xE11[iim1];
		    Np1b1[ii] = (z[iip1]+periodicYTranslation[2])*xE11[iip1];
		    Np2b1[ii] = (z[iip2]+periodicYTranslation[2])*xE11[iip2];
	
		    Nm2b3[ii] = (z[iim2]-periodicYTranslation[2])*xE13[iim2];
		    Nm1b3[ii] = (z[iim1]-periodicYTranslation[2])*xE13[iim1];
		    Np1b3[ii] = (z[iip1]+periodicYTranslation[2])*xE13[iip1];
		    Np2b3[ii] = (z[iip2]+periodicYTranslation[2])*xE13[iip2];
	
		    Nm2c1[ii] = (z[iim2]-periodicYTranslation[2])*xE21[iim2];
		    Nm1c1[ii] = (z[iim1]-periodicYTranslation[2])*xE21[iim1];
		    Np1c1[ii] = (z[iip1]+periodicYTranslation[2])*xE21[iip1];
		    Np2c1[ii] = (z[iip2]+periodicYTranslation[2])*xE21[iip2];
	
		    Nm2c3[ii] = (z[iim2]-periodicYTranslation[2])*xE23[iim2];
		    Nm1c3[ii] = (z[iim1]-periodicYTranslation[2])*xE23[iim1];
		    Np1c3[ii] = (z[iip1]+periodicYTranslation[2])*xE23[iip1];
		    Np2c3[ii] = (z[iip2]+periodicYTranslation[2])*xE23[iip2];
	

		}
	    }

	    derivY->calc1stDerivField_TPB(Va1, dVa12, Nm2a1, Nm1a1, Np1a1, Np2a1);
	    derivY->calc1stDerivField_TPB(Va3, dVa32, Nm2a3, Nm1a3, Np1a3, Np2a3);

	    derivY->calc1stDerivField_TPB(Vb1, dVb12, Nm2b1, Nm1b1, Np1b1, Np2b1);
	    derivY->calc1stDerivField_TPB(Vb3, dVb32, Nm2b3, Nm1b3, Np1b3, Np2b3);

	    derivY->calc1stDerivField_TPB(Vc1, dVc12, Nm2c1, Nm1c1, Np1c1, Np2c1);
	    derivY->calc1stDerivField_TPB(Vc3, dVc32, Nm2c3, Nm1c3, Np1c3, Np2c3);


	    delete[] Nm2a1;
	    delete[] Nm1a1;
	    delete[] Np1a1;
	    delete[] Np2a1;

	    delete[] Nm2a3;
	    delete[] Nm1a3;
	    delete[] Np1a3;
	    delete[] Np2a3;

	    delete[] Nm2b1;
	    delete[] Nm1b1;
	    delete[] Np1b1;
	    delete[] Np2b1;

	    delete[] Nm2b3;
	    delete[] Nm1b3;
	    delete[] Np1b3;
	    delete[] Np2b3;

	    delete[] Nm2c1;
	    delete[] Nm1c1;
	    delete[] Np1c1;
	    delete[] Np2c1;

	    delete[] Nm2c3;
	    delete[] Nm1c3;
	    delete[] Np1c3;
	    delete[] Np2c3;

	}else{
	    derivY->calc1stDerivField(Va1, dVa12);
	    derivY->calc1stDerivField(Va3, dVa32);
	    derivY->calc1stDerivField(Vb1, dVb12);
	    derivY->calc1stDerivField(Vb3, dVb32);
	    derivY->calc1stDerivField(Vc1, dVc12);
  	    derivY->calc1stDerivField(Vc3, dVc32);
	}

	//Start doing the E1 derivatives...
	double *tempX5, *tempX6, *tempX7, *tempX8;
	double *tempX9, *tempX10, *tempX11, *tempX12;
	c2d->allocX(tempX5); c2d->allocX(tempX9);
	c2d->allocX(tempX6); c2d->allocX(tempX10);
	c2d->allocX(tempX7); c2d->allocX(tempX11);
	c2d->allocX(tempX8); c2d->allocX(tempX12);

	double *y1, *z1;
	double *xE12_1, *xE13_1;
	double *xE22_1, *xE23_1;

	c2d->transposeY2X_MajorIndex(Va2, tempX1);
	c2d->transposeY2X_MajorIndex(Va3, tempX2);
	c2d->transposeY2X_MajorIndex(Vb2, tempX3);
	c2d->transposeY2X_MajorIndex(Vb3, tempX4);
	c2d->transposeY2X_MajorIndex(Vc2, tempX5);
	c2d->transposeY2X_MajorIndex(Vc3, tempX6);


	if(periodicX){
	    
	    c2d->allocX(y1);
	    c2d->allocX(z1);
	    c2d->transposeY2X_MajorIndex(y, y1);
	    c2d->transposeY2X_MajorIndex(z, z1);

	    c2d->allocX(xE12_1);
	    c2d->allocX(xE13_1);
	    c2d->allocX(xE22_1);
	    c2d->allocX(xE23_1);
	    c2d->transposeY2X_MajorIndex(xE12, xE12_1);
	    c2d->transposeY2X_MajorIndex(xE13, xE13_1);
	    c2d->transposeY2X_MajorIndex(xE22, xE22_1);
	    c2d->transposeY2X_MajorIndex(xE23, xE23_1);

	    double *Nm2a2, *Nm1a2, *Np1a2, *Np2a2;
            Nm2a2 = new double[pxSize[1]*pxSize[2]];
            Nm1a2 = new double[pxSize[1]*pxSize[2]];
            Np1a2 = new double[pxSize[1]*pxSize[2]];
            Np2a2 = new double[pxSize[1]*pxSize[2]];

	    double *Nm2a3, *Nm1a3, *Np1a3, *Np2a3;
            Nm2a3 = new double[pxSize[1]*pxSize[2]];
            Nm1a3 = new double[pxSize[1]*pxSize[2]];
            Np1a3 = new double[pxSize[1]*pxSize[2]];
            Np2a3 = new double[pxSize[1]*pxSize[2]];

	    double *Nm2b2, *Nm1b2, *Np1b2, *Np2b2;
	    Nm2b2 = new double[pxSize[1]*pxSize[2]];
	    Nm1b2 = new double[pxSize[1]*pxSize[2]];
	    Np1b2 = new double[pxSize[1]*pxSize[2]];
	    Np2b2 = new double[pxSize[1]*pxSize[2]];

	    double *Nm2b3, *Nm1b3, *Np1b3, *Np2b3;
	    Nm2b3 = new double[pxSize[1]*pxSize[2]];
	    Nm1b3 = new double[pxSize[1]*pxSize[2]];
	    Np1b3 = new double[pxSize[1]*pxSize[2]];
	    Np2b3 = new double[pxSize[1]*pxSize[2]];

	    double *Nm2c2, *Nm1c2, *Np1c2, *Np2c2;
	    Nm2c2 = new double[pxSize[1]*pxSize[2]];
	    Nm1c2 = new double[pxSize[1]*pxSize[2]];
	    Np1c2 = new double[pxSize[1]*pxSize[2]];
	    Np2c2 = new double[pxSize[1]*pxSize[2]];

	    double *Nm2c3, *Nm1c3, *Np1c3, *Np2c3;
	    Nm2c3 = new double[pxSize[1]*pxSize[2]];
	    Nm1c3 = new double[pxSize[1]*pxSize[2]];
	    Np1c3 = new double[pxSize[1]*pxSize[2]];
	    Np2c3 = new double[pxSize[1]*pxSize[2]];

            FOR_Z_XPEN{
                FOR_Y_XPEN{
                    int ii = k*pxSize[1] + j;

                    int iim2 = k*pxSize[0]*pxSize[1] + j*pxSize[0] + pxSize[0]-2;
                    int iim1 = k*pxSize[0]*pxSize[1] + j*pxSize[0] + pxSize[0]-1;
                    int iip1 = k*pxSize[0]*pxSize[1] + j*pxSize[0] + 0;
                    int iip2 = k*pxSize[0]*pxSize[1] + j*pxSize[0] + 1;

                    Nm2a2[ii] = (y1[iim2]-periodicXTranslation[1])*xE12_1[iim2];
                    Nm1a2[ii] = (y1[iim1]-periodicXTranslation[1])*xE12_1[iim1];
                    Np1a2[ii] = (y1[iip1]+periodicXTranslation[1])*xE12_1[iip1];
                    Np2a2[ii] = (y1[iip2]+periodicXTranslation[1])*xE12_1[iip2];
              
		    Nm2a3[ii] = (y1[iim2]-periodicXTranslation[1])*xE13_1[iim2];
                    Nm1a3[ii] = (y1[iim1]-periodicXTranslation[1])*xE13_1[iim1];
                    Np1a3[ii] = (y1[iip1]+periodicXTranslation[1])*xE13_1[iip1];
                    Np2a3[ii] = (y1[iip2]+periodicXTranslation[1])*xE13_1[iip2];
             
		    Nm2b2[ii] = (z1[iim2]-periodicXTranslation[2])*xE12_1[iim2];
                    Nm1b2[ii] = (z1[iim1]-periodicXTranslation[2])*xE12_1[iim1];
                    Np1b2[ii] = (z1[iip1]+periodicXTranslation[2])*xE12_1[iip1];
                    Np2b2[ii] = (z1[iip2]+periodicXTranslation[2])*xE12_1[iip2];
 
		    Nm2b3[ii] = (z1[iim2]-periodicXTranslation[2])*xE13_1[iim2];
                    Nm1b3[ii] = (z1[iim1]-periodicXTranslation[2])*xE13_1[iim1];
                    Np1b3[ii] = (z1[iip1]+periodicXTranslation[2])*xE13_1[iip1];
                    Np2b3[ii] = (z1[iip2]+periodicXTranslation[2])*xE13_1[iip2];
 
		    Nm2c2[ii] = (z1[iim2]-periodicXTranslation[2])*xE22_1[iim2];
                    Nm1c2[ii] = (z1[iim1]-periodicXTranslation[2])*xE22_1[iim1];
                    Np1c2[ii] = (z1[iip1]+periodicXTranslation[2])*xE22_1[iip1];
                    Np2c2[ii] = (z1[iip2]+periodicXTranslation[2])*xE22_1[iip2];

		    Nm2c3[ii] = (z1[iim2]-periodicXTranslation[2])*xE23_1[iim2];
                    Nm1c3[ii] = (z1[iim1]-periodicXTranslation[2])*xE23_1[iim1];
                    Np1c3[ii] = (z1[iip1]+periodicXTranslation[2])*xE23_1[iip1];
                    Np2c3[ii] = (z1[iip2]+periodicXTranslation[2])*xE23_1[iip2];
 
		}
            }

            derivX->calc1stDerivField_TPB(tempX1, tempX7,  Nm2a2, Nm1a2, Np1a2, Np2a2);
            derivX->calc1stDerivField_TPB(tempX2, tempX8,  Nm2a3, Nm1a3, Np1a3, Np2a3);
            derivX->calc1stDerivField_TPB(tempX3, tempX9,  Nm2b2, Nm1b2, Np1b2, Np2b2);
            derivX->calc1stDerivField_TPB(tempX4, tempX10, Nm2b3, Nm1b3, Np1b3, Np2b3);
            derivX->calc1stDerivField_TPB(tempX5, tempX11, Nm2c2, Nm1c2, Np1c2, Np2c2);
            derivX->calc1stDerivField_TPB(tempX6, tempX12, Nm2c3, Nm1c3, Np1c3, Np2c3);

	    c2d->deallocXYZ(y1);
	    c2d->deallocXYZ(z1);
	    c2d->deallocXYZ(xE12_1);
	    c2d->deallocXYZ(xE13_1);
	    c2d->deallocXYZ(xE22_1);
	    c2d->deallocXYZ(xE23_1);

            delete[] Nm2a2;
            delete[] Nm1a2;
            delete[] Np1a2;
            delete[] Np2a2;

            delete[] Nm2a3;
            delete[] Nm1a3;
            delete[] Np1a3;
            delete[] Np2a3;

            delete[] Nm2b2;
            delete[] Nm1b2;
            delete[] Np1b2;
            delete[] Np2b2;

            delete[] Nm2b3;
            delete[] Nm1b3;
            delete[] Np1b3;
            delete[] Np2b3;

            delete[] Nm2c2;
            delete[] Nm1c2;
            delete[] Np1c2;
            delete[] Np2c2;

            delete[] Nm2c3;
            delete[] Nm1c3;
            delete[] Np1c3;
            delete[] Np2c3;

	}else{
	    derivX->calc1stDerivField(tempX1, tempX7);
	    derivX->calc1stDerivField(tempX2, tempX8);
	    derivX->calc1stDerivField(tempX3, tempX9);
	    derivX->calc1stDerivField(tempX4, tempX10);
	    derivX->calc1stDerivField(tempX5, tempX11);
	    derivX->calc1stDerivField(tempX6, tempX12);
	}

	c2d->transposeX2Y_MajorIndex(tempX7, dVa21);
	c2d->transposeX2Y_MajorIndex(tempX8, dVa31);
	c2d->transposeX2Y_MajorIndex(tempX9, dVb21);
	c2d->transposeX2Y_MajorIndex(tempX10, dVb31);
	c2d->transposeX2Y_MajorIndex(tempX11, dVc21);
	c2d->transposeX2Y_MajorIndex(tempX12, dVc31);

	//Start doing the E3 derivatives...
	double *tempZ5, *tempZ6, *tempZ7, *tempZ8;
	double *tempZ9, *tempZ10, *tempZ11, *tempZ12;
	c2d->allocZ(tempZ5); c2d->allocZ(tempZ9);
	c2d->allocZ(tempZ6); c2d->allocZ(tempZ10);
	c2d->allocZ(tempZ7); c2d->allocZ(tempZ11);
	c2d->allocZ(tempZ8); c2d->allocZ(tempZ12);

	double *y3, *z3;
	double *xE11_3, *xE12_3;
	double *xE21_3, *xE22_3;



	c2d->transposeY2Z_MajorIndex(Va1, tempZ1);
	c2d->transposeY2Z_MajorIndex(Va2, tempZ2);
	c2d->transposeY2Z_MajorIndex(Vb1, tempZ3);
	c2d->transposeY2Z_MajorIndex(Vb2, tempZ4);
	c2d->transposeY2Z_MajorIndex(Vc1, tempZ5);
	c2d->transposeY2Z_MajorIndex(Vc2, tempZ6);

	if(periodicZ){
	
 	    c2d->allocZ(y3);
	    c2d->allocZ(z3);
	    c2d->transposeY2Z_MajorIndex(y, y3);
	    c2d->transposeY2Z_MajorIndex(z, z3);

	    c2d->allocZ(xE11_3);
	    c2d->allocZ(xE12_3);
	    c2d->allocZ(xE21_3);
	    c2d->allocZ(xE22_3);
	    c2d->transposeY2Z_MajorIndex(xE11, xE11_3);
	    c2d->transposeY2Z_MajorIndex(xE12, xE12_3);
	    c2d->transposeY2Z_MajorIndex(xE21, xE21_3);
	    c2d->transposeY2Z_MajorIndex(xE22, xE22_3);

	    double *Nm2a1, *Nm1a1, *Np1a1, *Np2a1;
            Nm2a1 = new double[pzSize[1]*pzSize[0]];
            Nm1a1 = new double[pzSize[1]*pzSize[0]];
            Np1a1 = new double[pzSize[1]*pzSize[0]];
            Np2a1 = new double[pzSize[1]*pzSize[0]];

	    double *Nm2a2, *Nm1a2, *Np1a2, *Np2a2;
            Nm2a2 = new double[pzSize[1]*pzSize[0]];
            Nm1a2 = new double[pzSize[1]*pzSize[0]];
            Np1a2 = new double[pzSize[1]*pzSize[0]];
            Np2a2 = new double[pzSize[1]*pzSize[0]];

	    double *Nm2b1, *Nm1b1, *Np1b1, *Np2b1;
            Nm2b1 = new double[pzSize[1]*pzSize[0]];
            Nm1b1 = new double[pzSize[1]*pzSize[0]];
            Np1b1 = new double[pzSize[1]*pzSize[0]];
            Np2b1 = new double[pzSize[1]*pzSize[0]];

	    double *Nm2b2, *Nm1b2, *Np1b2, *Np2b2;
            Nm2b2 = new double[pzSize[1]*pzSize[0]];
            Nm1b2 = new double[pzSize[1]*pzSize[0]];
            Np1b2 = new double[pzSize[1]*pzSize[0]];
            Np2b2 = new double[pzSize[1]*pzSize[0]];

	    double *Nm2c1, *Nm1c1, *Np1c1, *Np2c1;
            Nm2c1 = new double[pzSize[1]*pzSize[0]];
            Nm1c1 = new double[pzSize[1]*pzSize[0]];
            Np1c1 = new double[pzSize[1]*pzSize[0]];
            Np2c1 = new double[pzSize[1]*pzSize[0]];

	    double *Nm2c2, *Nm1c2, *Np1c2, *Np2c2;
            Nm2c2 = new double[pzSize[1]*pzSize[0]];
            Nm1c2 = new double[pzSize[1]*pzSize[0]];
            Np1c2 = new double[pzSize[1]*pzSize[0]];
            Np2c2 = new double[pzSize[1]*pzSize[0]];

            FOR_Y_ZPEN{
                FOR_X_ZPEN{
                    int ii = j*pzSize[0] + i;

                    int iim2 = j*pzSize[0]*pzSize[2] + i*pzSize[2] + pzSize[2]-2;
                    int iim1 = j*pzSize[0]*pzSize[2] + i*pzSize[2] + pzSize[2]-1;
                    int iip1 = j*pzSize[0]*pzSize[2] + i*pzSize[2] + 0;
                    int iip2 = j*pzSize[0]*pzSize[2] + i*pzSize[2] + 1;

                    Nm2a1[ii] = (y3[iim2]-periodicZTranslation[1])*xE11_3[iim2];
                    Nm1a1[ii] = (y3[iim1]-periodicZTranslation[1])*xE11_3[iim1];
                    Np1a1[ii] = (y3[iip1]+periodicZTranslation[1])*xE11_3[iip1];
                    Np2a1[ii] = (y3[iip2]+periodicZTranslation[1])*xE11_3[iip2];

                    Nm2a2[ii] = (y3[iim2]-periodicZTranslation[1])*xE12_3[iim2];
                    Nm1a2[ii] = (y3[iim1]-periodicZTranslation[1])*xE12_3[iim1];
                    Np1a2[ii] = (y3[iip1]+periodicZTranslation[1])*xE12_3[iip1];
                    Np2a2[ii] = (y3[iip2]+periodicZTranslation[1])*xE12_3[iip2];

                    Nm2b1[ii] = (z3[iim2]-periodicZTranslation[2])*xE11_3[iim2];
                    Nm1b1[ii] = (z3[iim1]-periodicZTranslation[2])*xE11_3[iim1];
                    Np1b1[ii] = (z3[iip1]+periodicZTranslation[2])*xE11_3[iip1];
                    Np2b1[ii] = (z3[iip2]+periodicZTranslation[2])*xE11_3[iip2];

                    Nm2b2[ii] = (z3[iim2]-periodicZTranslation[2])*xE12_3[iim2];
                    Nm1b2[ii] = (z3[iim1]-periodicZTranslation[2])*xE12_3[iim1];
                    Np1b2[ii] = (z3[iip1]+periodicZTranslation[2])*xE12_3[iip1];
                    Np2b2[ii] = (z3[iip2]+periodicZTranslation[2])*xE12_3[iip2];

                    Nm2c1[ii] = (z3[iim2]-periodicZTranslation[2])*xE21_3[iim2];
                    Nm1c1[ii] = (z3[iim1]-periodicZTranslation[2])*xE21_3[iim1];
                    Np1c1[ii] = (z3[iip1]+periodicZTranslation[2])*xE21_3[iip1];
                    Np2c1[ii] = (z3[iip2]+periodicZTranslation[2])*xE21_3[iip2];

                    Nm2c2[ii] = (z3[iim2]-periodicZTranslation[2])*xE22_3[iim2];
                    Nm1c2[ii] = (z3[iim1]-periodicZTranslation[2])*xE22_3[iim1];
                    Np1c2[ii] = (z3[iip1]+periodicZTranslation[2])*xE22_3[iip1];
                    Np2c2[ii] = (z3[iip2]+periodicZTranslation[2])*xE22_3[iip2];

                }
            }

            derivZ->calc1stDerivField_TPB(tempZ1, tempZ7,  Nm2a1, Nm1a1, Np1a1, Np2a1);
            derivZ->calc1stDerivField_TPB(tempZ2, tempZ8,  Nm2a2, Nm1a2, Np1a2, Np2a2);
            derivZ->calc1stDerivField_TPB(tempZ3, tempZ9,  Nm2b1, Nm1b1, Np1b1, Np2b1);
            derivZ->calc1stDerivField_TPB(tempZ4, tempZ10, Nm2b2, Nm1b2, Np1b2, Np2b2);
            derivZ->calc1stDerivField_TPB(tempZ5, tempZ11, Nm2c1, Nm1c1, Np1c1, Np2c1);
	    derivZ->calc1stDerivField_TPB(tempZ6, tempZ12, Nm2c2, Nm1c2, Np1c2, Np2c2);


	    c2d->deallocXYZ(y3);
	    c2d->deallocXYZ(z3);
	    c2d->deallocXYZ(xE11_3);
	    c2d->deallocXYZ(xE12_3);
	    c2d->deallocXYZ(xE21_3);
	    c2d->deallocXYZ(xE22_3);

            delete[] Nm2a1;
            delete[] Nm1a1;
            delete[] Np1a1;
            delete[] Np2a1;

            delete[] Nm2a2;
            delete[] Nm1a2;
            delete[] Np1a2;
            delete[] Np2a2;

            delete[] Nm2b1;
            delete[] Nm1b1;
            delete[] Np1b1;
            delete[] Np2b1;

            delete[] Nm2b2;
            delete[] Nm1b2;
            delete[] Np1b2;
            delete[] Np2b2;

            delete[] Nm2c1;
            delete[] Nm1c1;
            delete[] Np1c1;
            delete[] Np2c1;

            delete[] Nm2c2;
            delete[] Nm1c2;
            delete[] Np1c2;
            delete[] Np2c2;

	}else{
	    derivZ->calc1stDerivField(tempZ1, tempZ7);
	    derivZ->calc1stDerivField(tempZ2, tempZ8);
	    derivZ->calc1stDerivField(tempZ3, tempZ9);
	    derivZ->calc1stDerivField(tempZ4, tempZ10);
	    derivZ->calc1stDerivField(tempZ5, tempZ11);
	    derivZ->calc1stDerivField(tempZ6, tempZ12);
	}

	c2d->transposeZ2Y_MajorIndex(tempZ7, dVa13);
	c2d->transposeZ2Y_MajorIndex(tempZ8, dVa23);
	c2d->transposeZ2Y_MajorIndex(tempZ9, dVb13);
	c2d->transposeZ2Y_MajorIndex(tempZ10, dVb23);
	c2d->transposeZ2Y_MajorIndex(tempZ11, dVc13);
	c2d->transposeZ2Y_MajorIndex(tempZ12, dVc23);


	c2d->deallocXYZ(xE11);
	c2d->deallocXYZ(xE12);
	c2d->deallocXYZ(xE13);
	c2d->deallocXYZ(xE21);
	c2d->deallocXYZ(xE22);
	c2d->deallocXYZ(xE23);

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
	    J33[ip] = dVa12[ip] - dVa21[ip];
	}

	//Free up all of this space here...
	
	c2d->deallocXYZ(Va1);
	c2d->deallocXYZ(Va2);
	c2d->deallocXYZ(Va3);
	c2d->deallocXYZ(Vb1);
	c2d->deallocXYZ(Vb2);
	c2d->deallocXYZ(Vb3);
	c2d->deallocXYZ(Vc1);
	c2d->deallocXYZ(Vc2);
	c2d->deallocXYZ(Vc3);

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
	if(periodicY){
	    double *Nm2, *Nm1, *Np1, *Np2;
	    Nm2 = new double[pySize[0]*pySize[2]];
	    Nm1 = new double[pySize[0]*pySize[2]];
	    Np1 = new double[pySize[0]*pySize[2]];
	    Np2 = new double[pySize[0]*pySize[2]];

	    FOR_X_YPEN{
		FOR_Z_YPEN{
		    int ii = i*pySize[2] + k;

		    int iim2 = i*pySize[2]*pySize[1] + k*pySize[1] + pySize[1]-2;			     
		    int iim1 = i*pySize[2]*pySize[1] + k*pySize[1] + pySize[1]-1;		
		    int iip1 = i*pySize[2]*pySize[1] + k*pySize[1] + 0;		
		    int iip2 = i*pySize[2]*pySize[1] + k*pySize[1] + 1;		
	
		    Nm2[ii] = (x[iim2]-periodicYTranslation[0])*J21[iim2] + 
			      (y[iim2]-periodicYTranslation[1])*J22[iim2] +
			      (z[iim2]-periodicYTranslation[2])*J23[iim2];

		    Nm1[ii] = (x[iim1]-periodicYTranslation[0])*J21[iim1] + 
			      (y[iim1]-periodicYTranslation[1])*J22[iim1] +
			      (z[iim1]-periodicYTranslation[2])*J23[iim1];

		    Np1[ii] = (x[iip1]+periodicYTranslation[0])*J21[iip1] + 
			      (y[iip1]+periodicYTranslation[1])*J22[iip1] +
			      (z[iip1]+periodicYTranslation[2])*J23[iip1];

		    Np2[ii] = (x[iip2]+periodicYTranslation[0])*J21[iip2] + 
			      (y[iip2]+periodicYTranslation[1])*J22[iip2] +
			      (z[iip2]+periodicYTranslation[2])*J23[iip2];

		}
	    }

	    derivY->calc1stDerivField_TPB(preJdet2, Jdet2, Nm2, Nm1, Np1, Np2);

	    delete[] Nm2;
	    delete[] Nm1;
	    delete[] Np1;
	    delete[] Np2;
	
	}else{
	    derivY->calc1stDerivField(preJdet2, Jdet2);
	}

	//Compute the E1 component....
	c2d->transposeY2X_MajorIndex(preJdet1, tempX1);

	if(periodicX){

	    c2d->transposeY2X_MajorIndex(x, tempX2);
	    c2d->transposeY2X_MajorIndex(y, tempX3);
	    c2d->transposeY2X_MajorIndex(z, tempX4);

	    c2d->transposeY2X_MajorIndex(J11, tempX6);
	    c2d->transposeY2X_MajorIndex(J12, tempX7);
	    c2d->transposeY2X_MajorIndex(J13, tempX8);

	    double *Nm2, *Nm1, *Np1, *Np2;
	    Nm2 = new double[pxSize[1]*pxSize[2]];
	    Nm1 = new double[pxSize[1]*pxSize[2]];
	    Np1 = new double[pxSize[1]*pxSize[2]];
	    Np2 = new double[pxSize[1]*pxSize[2]];

	    FOR_Z_XPEN{
                FOR_Y_XPEN{
                    int ii = k*pxSize[1] + j;

                    int iim2 = k*pxSize[0]*pxSize[1] + j*pxSize[0] + pxSize[0]-2;
                    int iim1 = k*pxSize[0]*pxSize[1] + j*pxSize[0] + pxSize[0]-1;
                    int iip1 = k*pxSize[0]*pxSize[1] + j*pxSize[0] + 0;
                    int iip2 = k*pxSize[0]*pxSize[1] + j*pxSize[0] + 1;

		    Nm2[ii] = (tempX2[iim2]-periodicXTranslation[0])*tempX6[iim2] + 
			      (tempX3[iim2]-periodicXTranslation[1])*tempX7[iim2] +
			      (tempX4[iim2]-periodicXTranslation[2])*tempX8[iim2];

		    Nm1[ii] = (tempX2[iim1]-periodicXTranslation[0])*tempX6[iim1] + 
			      (tempX3[iim1]-periodicXTranslation[1])*tempX7[iim1] +
			      (tempX4[iim1]-periodicXTranslation[2])*tempX8[iim1];

		    Np1[ii] = (tempX2[iip1]+periodicXTranslation[0])*tempX6[iip1] + 
			      (tempX3[iip1]+periodicXTranslation[1])*tempX7[iip1] +
			      (tempX4[iip1]+periodicXTranslation[2])*tempX8[iip1];

		    Np2[ii] = (tempX2[iip2]+periodicXTranslation[0])*tempX6[iip2] + 
			      (tempX3[iip2]+periodicXTranslation[1])*tempX7[iip2] +
			      (tempX4[iip2]+periodicXTranslation[2])*tempX8[iip2];

		}
	    }	


	    derivX->calc1stDerivField_TPB(tempX1, tempX5, Nm2, Nm1, Np1, Np2);

	    delete[] Nm2;
	    delete[] Nm1;
	    delete[] Np1;
	    delete[] Np2;

	}else{
	    derivX->calc1stDerivField(tempX1, tempX5);
	}

	c2d->transposeX2Y_MajorIndex(tempX5, Jdet1);

 	//Compute the E3 component....
	c2d->transposeY2Z_MajorIndex(preJdet3, tempZ1);

	if(periodicZ){

	    c2d->transposeY2Z_MajorIndex(x, tempZ2);
	    c2d->transposeY2Z_MajorIndex(y, tempZ3);
	    c2d->transposeY2Z_MajorIndex(z, tempZ4);

	    c2d->transposeY2Z_MajorIndex(J31, tempZ6);
	    c2d->transposeY2Z_MajorIndex(J32, tempZ7);
	    c2d->transposeY2Z_MajorIndex(J33, tempZ8);

	    double *Nm2, *Nm1, *Np1, *Np2;
	    Nm2 = new double[pzSize[0]*pzSize[1]];
	    Nm1 = new double[pzSize[0]*pzSize[1]];
	    Np1 = new double[pzSize[0]*pzSize[1]];
	    Np2 = new double[pzSize[0]*pzSize[1]];

            FOR_Y_ZPEN{
                FOR_X_ZPEN{
                    int ii = j*pzSize[0] + i;

                    int iim2 = j*pzSize[0]*pzSize[2] + i*pzSize[2] + pzSize[2]-2;
                    int iim1 = j*pzSize[0]*pzSize[2] + i*pzSize[2] + pzSize[2]-1;
                    int iip1 = j*pzSize[0]*pzSize[2] + i*pzSize[2] + 0;
                    int iip2 = j*pzSize[0]*pzSize[2] + i*pzSize[2] + 1;

		    Nm2[ii] = (tempZ2[iim2]-periodicZTranslation[0])*tempZ6[iim2] + 
			      (tempZ3[iim2]-periodicZTranslation[1])*tempZ7[iim2] +
			      (tempZ4[iim2]-periodicZTranslation[2])*tempZ8[iim2];

		    Nm1[ii] = (tempZ2[iim1]-periodicZTranslation[0])*tempZ6[iim1] + 
			      (tempZ3[iim1]-periodicZTranslation[1])*tempZ7[iim1] +
			      (tempZ4[iim1]-periodicZTranslation[2])*tempZ8[iim1];

		    Np1[ii] = (tempZ2[iip1]+periodicZTranslation[0])*tempZ6[iip1] + 
			      (tempZ3[iip1]+periodicZTranslation[1])*tempZ7[iip1] +
			      (tempZ4[iip1]+periodicZTranslation[2])*tempZ8[iip1];

		    Np2[ii] = (tempZ2[iip2]+periodicZTranslation[0])*tempZ6[iip2] + 
			      (tempZ3[iip2]+periodicZTranslation[1])*tempZ7[iip2] +
			      (tempZ4[iip2]+periodicZTranslation[2])*tempZ8[iip2];

		}
	    }

	    derivZ->calc1stDerivField_TPB(tempZ1, tempZ5, Nm2, Nm1, Np1, Np2);

	    delete[] Nm2;
	    delete[] Nm1;
	    delete[] Np1;
	    delete[] Np2;

	}else{
	    derivZ->calc1stDerivField(tempZ1, tempZ5);
	}

	c2d->transposeZ2Y_MajorIndex(tempZ5, Jdet3);

	//Compute the Jacobian determinant...
	FOR_XYZ_YPEN{
	    J[ip] = 1/((1.0/3.0)*(Jdet1[ip] + Jdet2[ip] + Jdet3[ip]));
	    J11[ip] = J[ip]*J11[ip];
	    J12[ip] = J[ip]*J12[ip];
	    J13[ip] = J[ip]*J13[ip];
	    J21[ip] = J[ip]*J21[ip];
	    J22[ip] = J[ip]*J22[ip];
	    J23[ip] = J[ip]*J23[ip];
	    J31[ip] = J[ip]*J31[ip];
	    J32[ip] = J[ip]*J32[ip];
	    J33[ip] = J[ip]*J33[ip];
	}

	IF_RANK0 cout << "check" << endl;

	//Copy data over to solver to make simpler to access
	cs->J = J;
	cs->J11 = J11;
	cs->J12 = J12;
	cs->J13 = J13;
	cs->J21 = J21;
	cs->J22 = J22;
	cs->J23 = J23;
	cs->J31 = J31;
	cs->J32 = J32;
	cs->J33 = J33;

	IF_RANK0 cout << "done!" << endl;
	
        getRange(J, "J", pySize[0], pySize[1], pySize[2], mpiRank);
        getRange(J11, "J11", pySize[0], pySize[1], pySize[2], mpiRank);
        getRange(J12, "J12", pySize[0], pySize[1], pySize[2], mpiRank);
        getRange(J13, "J13", pySize[0], pySize[1], pySize[2], mpiRank);
        getRange(J21, "J21", pySize[0], pySize[1], pySize[2], mpiRank);
        getRange(J22, "J22", pySize[0], pySize[1], pySize[2], mpiRank);
        getRange(J23, "J23", pySize[0], pySize[1], pySize[2], mpiRank);
        getRange(J31, "J31", pySize[0], pySize[1], pySize[2], mpiRank);
        getRange(J32, "J32", pySize[0], pySize[1], pySize[2], mpiRank);
        getRange(J33, "J33", pySize[0], pySize[1], pySize[2], mpiRank);

	

	IF_RANK0 cout << " > Checking the values of the metric indentities, values should be small" << endl;

	double *I1_1, *I1_2, *I1_3;
	double *I2_1, *I2_2, *I2_3;
	double *I3_1, *I3_2, *I3_3;

	c2d->allocY(I1_1);
	c2d->allocY(I1_2);
	c2d->allocY(I1_3);
	c2d->allocY(I2_1);
	c2d->allocY(I2_2);
	c2d->allocY(I2_3);
	c2d->allocY(I3_1);
	c2d->allocY(I3_2);
	c2d->allocY(I3_3);

	derivY->calc1stDerivField(J21, I1_2);
	derivY->calc1stDerivField(J22, I2_2);
	derivY->calc1stDerivField(J23, I3_2);


        c2d->transposeY2X_MajorIndex(J11, tempX1);
        c2d->transposeY2X_MajorIndex(J12, tempX2);
        c2d->transposeY2X_MajorIndex(J13, tempX3);

	derivX->calc1stDerivField(tempX1, tempX4);
	derivX->calc1stDerivField(tempX2, tempX5);
	derivX->calc1stDerivField(tempX3, tempX6);

	c2d->transposeX2Y_MajorIndex(tempX4, I1_1);
	c2d->transposeX2Y_MajorIndex(tempX5, I2_1);
	c2d->transposeX2Y_MajorIndex(tempX6, I3_1);


        c2d->transposeY2Z_MajorIndex(J31, tempZ1);
        c2d->transposeY2Z_MajorIndex(J32, tempZ2);
        c2d->transposeY2Z_MajorIndex(J33, tempZ3);

	derivZ->calc1stDerivField(tempZ1, tempZ4);
	derivZ->calc1stDerivField(tempZ2, tempZ5);
	derivZ->calc1stDerivField(tempZ3, tempZ6);

	c2d->transposeZ2Y_MajorIndex(tempZ4, I1_3);
	c2d->transposeZ2Y_MajorIndex(tempZ5, I2_3);
	c2d->transposeZ2Y_MajorIndex(tempZ6, I3_3);

	FOR_XYZ_YPEN{
	    I1_1[ip] = I1_1[ip] + I1_2[ip] + I1_3[ip];
	    I2_1[ip] = I2_1[ip] + I2_2[ip] + I2_3[ip];
	    I3_1[ip] = I3_1[ip] + I3_2[ip] + I3_3[ip];
	}

	getRange(I1_1, "Metric Identity 1", pySize[0], pySize[1], pySize[2], mpiRank);
	getRange(I2_1, "Metric Identity 2", pySize[0], pySize[1], pySize[2], mpiRank);
	getRange(I3_1, "Metric Identity 3", pySize[0], pySize[1], pySize[2], mpiRank);

	c2d->deallocXYZ(I1_1);
	c2d->deallocXYZ(I2_1);
	c2d->deallocXYZ(I3_1);
	c2d->deallocXYZ(I1_2);
	c2d->deallocXYZ(I2_2);
	c2d->deallocXYZ(I3_2);
	c2d->deallocXYZ(I1_3);
	c2d->deallocXYZ(I2_3);
	c2d->deallocXYZ(I3_3);

	//Free up all of the spaces we've been using...
	c2d->deallocXYZ(tempX1);
	c2d->deallocXYZ(tempX2);
	c2d->deallocXYZ(tempX3);
	c2d->deallocXYZ(tempX4);
	c2d->deallocXYZ(tempX5);
	c2d->deallocXYZ(tempX6);
	c2d->deallocXYZ(tempX7);
	c2d->deallocXYZ(tempX8);
	c2d->deallocXYZ(tempX9);
	c2d->deallocXYZ(tempX10);
	c2d->deallocXYZ(tempX11);
	c2d->deallocXYZ(tempX12);
	c2d->deallocXYZ(tempZ1);
	c2d->deallocXYZ(tempZ2);
	c2d->deallocXYZ(tempZ3);
	c2d->deallocXYZ(tempZ4);
	c2d->deallocXYZ(tempZ5);
	c2d->deallocXYZ(tempZ6);
	c2d->deallocXYZ(tempZ7);
	c2d->deallocXYZ(tempZ8);
	c2d->deallocXYZ(tempZ9);
	c2d->deallocXYZ(tempZ10);
	c2d->deallocXYZ(tempZ11);
	c2d->deallocXYZ(tempZ12);

	c2d->deallocXYZ(preJdet1);
	c2d->deallocXYZ(preJdet2);
	c2d->deallocXYZ(preJdet3);

	c2d->deallocXYZ(Jdet1);
	c2d->deallocXYZ(Jdet2);
	c2d->deallocXYZ(Jdet3);

}

#endif
