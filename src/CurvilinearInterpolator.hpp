#ifndef _CURVILINEARINTEROLATORH_
#define _CURVILINEARINTEROLATORH_

#include "AbstractCSolver.hpp"
#include "AbstractSingleBlockMesh.hpp"
#include "Macros.hpp"
#include "Utils.hpp"
#include "Domain.hpp"
#include "BC.hpp"
#include <vector>

class CurvilinearInterpolator{

  public:

    int mpiRank;
    int pxSize[3], pySize[3], pzSize[3];
    int pxStart[3], pyStart[3], pzStart[3];
    int pxEnd[3], pyEnd[3], pzEnd[3];
    int Nx, Ny, Nz;
    bool transPeriodicX, transPeriodicY, transPeriodicZ;

    AbstractCSolver *cs;
    double (*pointList)[3];
    vector<int> icvList;
    vector<int> pointIndex;
    vector<double> pointListX;
    vector<double> pointListY;
    vector<double> pointListZ;
    int localPointFoundCount, globalPointFoundCount;

    //Weights for the components of the box to interpolate
    double (*Ni)[8];

    ~CurvilinearInterpolator(){
	delete[] Ni;
	icvList.clear();
	pointIndex.clear();
	pointListX.clear();
	pointListY.clear();
	pointListZ.clear();
    }

    CurvilinearInterpolator(AbstractCSolver *cs, double (*pointList)[3], int Npoints){

	//Getting geometry and solver data from the input objects
	this->cs = cs;
	this->pointList = pointList;
	this->Nx = cs->msh->Nx;
	this->Ny = cs->msh->Ny;
	this->Nz = cs->msh->Nz;

	this->transPeriodicX = cs->msh->transPeriodicX;
	this->transPeriodicY = cs->msh->transPeriodicY;
	this->transPeriodicZ = cs->msh->transPeriodicZ;

	cs->dom->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);
	mpiRank = cs->mpiRank;

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
		pointIndex.push_back(ip);
		pointListX.push_back(p[0]);
		pointListY.push_back(p[1]);
		pointListZ.push_back(p[2]);
	    }
	}

	MPI_Allreduce(&localPointFoundCount, &globalPointFoundCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 

        if(cs->mpiRank == 0){
	    cout << " > Found " << globalPointFoundCount << " of " << Npoints << " points for interpolation in the domain." << endl;
	}

	//allocate the local weights
	Ni = new double[localPointFoundCount][8];

	//For each of local found points
	for(int ii = 0; ii < localPointFoundCount; ii++){

		
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


	    //Just implement inverse distance weighting for now since we're just using it for visualization
	    
	    //Get the distance from the point to the corners of the box
	    double d2[8], p[3];
	    p[0] = pointListX[ii];
	    p[1] = pointListY[ii];
	    p[2] = pointListZ[ii];
	    for(int i = 0; i < 8; i++){
		d2[i] = 0.0;
		FOR_J3{
		    d2[i] += (box_p[i][j]-p[j])*(box_p[i][j]-p[j]); 
		}
	    }
	   /* INVERSE DISTANCE WEIGHTING...
	    double total_weight = 0.0;
	    for(int i = 0; i < 8; i++){
		total_weight += 1.0/d2[i];
	    }	       
 
	    for(int i = 0; i < 8; i++){
		Ni[ii][i] = (1.0/d2[i])/total_weight;
	    }
	   */ 

	    // TRILINEAR INTERPOLATION
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


	    double f[8][3];

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
		    f[i][j] *= (1.0/8.0);
		}
	    }

	    bool done = false;
	    int iter = 0;
	    
	    double e1 = 0.0, e2 = 0.0, e3 = 0.0; 
	    while(!done){

		iter++;

	        double error_tol = 1E-12;

		double alpha = f[0][0] + f[1][0]*e1 + f[2][0]*e2 + f[3][0]*e3 + f[4][0]*e1*e2 + \
			       f[5][0]*e2*e3 + f[6][0]*e1*e3 + f[7][0]*e1*e2*e3 - pointListX[ii]; 
		double beta  = f[0][1] + f[1][1]*e1 + f[2][1]*e2 + f[3][1]*e3 + f[4][1]*e1*e2 + \
			       f[5][1]*e2*e3 + f[6][1]*e1*e3 + f[7][1]*e1*e2*e3 - pointListY[ii]; 
		double gamma = f[0][2] + f[1][2]*e1 + f[2][2]*e2 + f[3][2]*e3 + f[4][2]*e1*e2 + \
			       f[5][2]*e2*e3 + f[6][2]*e1*e3 + f[7][2]*e1*e2*e3 - pointListZ[ii]; 


	   	double J[3][3];

		J[0][0] = f[1][0] + f[4][0]*e2 + f[6][0]*e3 + f[7][0]*e2*e3;
		J[1][0] = f[1][1] + f[4][1]*e2 + f[6][1]*e3 + f[7][1]*e2*e3;
		J[2][0] = f[1][2] + f[4][2]*e2 + f[6][2]*e3 + f[7][2]*e2*e3;


		J[0][1] = f[2][0] + f[4][0]*e1 + f[5][0]*e3 + f[7][0]*e1*e3;
		J[1][1] = f[2][1] + f[4][1]*e1 + f[5][1]*e3 + f[7][1]*e1*e3;
		J[2][1] = f[2][2] + f[4][2]*e1 + f[5][2]*e3 + f[7][2]*e1*e3;


		J[0][2] = f[3][0] + f[5][0]*e2 + f[6][0]*e1 + f[7][0]*e1*e2;
		J[1][2] = f[3][1] + f[5][1]*e2 + f[6][1]*e1 + f[7][1]*e1*e2;
		J[2][2] = f[3][2] + f[5][2]*e2 + f[6][2]*e1 + f[7][2]*e1*e2;

		double Jinv[3][3];
		
		double A, B, C, D, E, F, G, H, I;
		A =  (J[1][1]*J[2][2] - J[1][2]*J[2][1]);
		B = -(J[1][0]*J[2][2] - J[1][2]*J[2][0]);
		C =  (J[1][0]*J[2][1] - J[1][1]*J[2][0]);
		D = -(J[0][1]*J[2][2] - J[0][2]*J[2][1]);
		E =  (J[0][0]*J[2][2] - J[0][2]*J[2][0]);
		F = -(J[0][0]*J[2][1] - J[0][1]*J[2][0]); 
		G =  (J[0][1]*J[1][2] - J[0][2]*J[1][1]);
		H = -(J[0][0]*J[1][2] - J[0][2]*J[1][0]);
		I =  (J[0][0]*J[1][1] - J[0][1]*J[1][0]);

		double Jdet = J[0][0]*A + J[1][1]*B + J[2][2]*C;

		
		Jinv[0][0] = A/Jdet;
		Jinv[1][0] = B/Jdet;
		Jinv[2][0] = C/Jdet;
		Jinv[0][1] = D/Jdet;
		Jinv[1][1] = E/Jdet;
		Jinv[2][1] = F/Jdet;
		Jinv[0][2] = G/Jdet;
		Jinv[1][2] = H/Jdet;
		Jinv[2][2] = I/Jdet;

		double delta[3];
		FOR_I3	delta[i] = alpha*Jinv[i][0] + beta*Jinv[i][1] + gamma*Jinv[i][2];
			
		double e1_new = e1 - delta[0];
		double e2_new = e2 - delta[1];
		double e3_new = e3 - delta[2];
 
		double eps[3];
		eps[0] = fabs(e1_new-e1); 
		eps[1] = fabs(e2_new-e2); 
		eps[2] = fabs(e3_new-e3); 

		if(eps[0] < error_tol && eps[1] < error_tol && eps[2] < error_tol){
		    done = true;
		}

		e1 = e1_new;
		e2 = e2_new;
		e3 = e3_new;

	    }

	    //Now we have to get and reorganize the weights to be in the right place
	    Ni[ii][7] = (1.0/8.0)*(1.0+e1)*(1.0+e2)*(1.0+e3);
	    Ni[ii][3] = (1.0/8.0)*(1.0-e1)*(1.0+e2)*(1.0+e3);
	    Ni[ii][2] = (1.0/8.0)*(1.0-e1)*(1.0+e2)*(1.0-e3);
	    Ni[ii][6] = (1.0/8.0)*(1.0+e1)*(1.0+e2)*(1.0-e3);
	    Ni[ii][5] = (1.0/8.0)*(1.0+e1)*(1.0-e2)*(1.0+e3);
	    Ni[ii][1] = (1.0/8.0)*(1.0-e1)*(1.0-e2)*(1.0+e3);
	    Ni[ii][0] = (1.0/8.0)*(1.0-e1)*(1.0-e2)*(1.0-e3);
	    Ni[ii][4] = (1.0/8.0)*(1.0+e1)*(1.0-e2)*(1.0-e3);

	}



	delete[] xHalo;
	delete[] yHalo;
	delete[] zHalo;

    }

    void getDataHalo(double *dataIn, double *&dataHalo);

    void getOrderedBlockData(int ip, int jp, int kp, double *dataHalo, double box_data[8]);

    void interpolateData(double *dataIn, double *interpedDataOut);  

};



#endif
