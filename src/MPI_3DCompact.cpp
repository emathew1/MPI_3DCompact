#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>

using namespace std;

#include "C2Decomp.hpp"
#include "Macros.hpp"
#include "Options.hpp"
#include "TimeStepping.hpp"
#include "Domain.hpp"
#include "BC.hpp" 

#include "AbstractCSolver.hpp"
#include "CurvilinearCSolver.hpp"

#include "AbstractSingleBlockMesh.hpp"
#include "AlgebraicSingleBlockMesh.hpp"

#include "AbstractRK.hpp"
#include "TVDRK3.hpp"

#include "CurvilinearInterpolator.hpp"

int main(int argc, char *argv[]){
   int ierr, totRank, mpiRank;

    //Initialize MPI
    ierr = MPI_Init( &argc, &argv);

    //Get the number of processes
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &totRank);

    //Get the local rank
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    IF_RANK0{
        cout << endl;
        cout << " --------------------" << endl;
    	cout << "  Starting up Solver " << endl;
    	cout << " --------------------" << endl;
	cout << endl;
	cout << " > Running on " << totRank << " cores!" << endl;
	cout << endl;
    }


    //Have the root rank parse the input file and broadcast it out...
    Options *opt = new Options(mpiRank);



    ////////////////////////////////////
    //Time Stepping info intialization//
    ////////////////////////////////////
    TimeStepping *ts = new TimeStepping(opt);


    ///////////////////////////
    //Boundary Condition Info//
    ///////////////////////////

    bool periodicBC[3];
    BC *bc = new BC(opt, periodicBC);


    //Get the dimesions of the grid so we can initialize the C2Decomp object...

    MPI_Offset disp;
    int TS, Nx, Ny, Nz;
    double time;
    if(opt->fromRestart || opt->onlyGridFromRestart){

	//If we need to read from a file, pull the dimensions from
	//the leading three doubles from the file...

	double cN[5];
	IF_RANK0{
	    FILE *ptr;
	    ptr = fopen(opt->filename.c_str(), "rb");
	    if(ptr == NULL){
		cout << "ERROR: Couldn't open file " << opt->filename << endl;
		MPI_Abort(MPI_COMM_WORLD, -10);
	    }else{
	        fread(cN, sizeof(double), 5, ptr);
	    }
	    fclose(ptr);
	}

	MPI_Bcast(cN, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	TS = (int)cN[0];
	time = cN[1];
	Nx = (int)cN[2];
	Ny = (int)cN[3];
	Nz = (int)cN[4];

	if(opt->fromRestart){
	    opt->timeStep = TS;
	    opt->time = time;
	}else{
	    opt->timeStep = 0;
	    opt->time = 0.0;
	}
	opt->Nx = Nx;
	opt->Ny = Ny;
	opt->Nz = Nz;

    }else{
	TS = 0;
	time = 0.0;
	opt->timeStep = TS;
	opt->time = time;
	Nx = opt->Nx;
	Ny = opt->Ny;
	Nz = opt->Nz;
	disp = 0;
    }

    /////////////////////////
    //Initialize the Domain//
    /////////////////////////
    //For curvilinear coordinates these should all correspond to the max xi, eta, and zeta values
    double Lx = 1.0, Ly = 1.0, Lz = 1.0;
    Domain *d = new Domain(opt, Lx, Ly, Lz);

    /////////////////////////////
    //Initializing Pencil Decomp//
    //////////////////////////////
 
    IF_RANK0 cout << endl << " > Initializing the pencil decomposition... " << endl;
    C2Decomp *c2d = new C2Decomp(opt->Nx, opt->Ny, opt->Nz, opt->pRow, opt->pCol, periodicBC);
    IF_RANK0 cout << endl << " > Handing some decomp info back to the domain object... " << endl;
    d->setPencilDecompInfo(c2d);

    //This is info needed to use the macros...
    int pxSize[3], pySize[3], pzSize[3];
    int pxStart[3], pyStart[3], pzStart[3];
    int pxEnd[3], pyEnd[3], pzEnd[3];
    d->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);


    /////////////////////////
    //Initialize the Solver//
    /////////////////////////
    AbstractCSolver *cs= new CurvilinearCSolver(c2d, d, bc, ts, opt);

    //Attach the mesh object to the solver...
    cs->msh = new AlgebraicSingleBlockMesh(c2d, cs, d, mpiRank);

    ///////////////////////////////////////////
    //Initialize Execution Loop and RK Method//
    ///////////////////////////////////////////
    AbstractRK *rk = new TVDRK3(cs);

    ///////////////////////////////
    //Set flow initial conditions//
    ///////////////////////////////

    bool fromRestart = opt->fromRestart;;

    if(!fromRestart){
        FOR_Z_YPEN{
            FOR_Y_YPEN{
                FOR_X_YPEN{

                    int ip = GETMAJIND_YPEN;

		    int ii = GETGLOBALXIND_YPEN;		
		    int jj = GETGLOBALYIND_YPEN;		
		    int kk = GETGLOBALZIND_YPEN;

		    double x  = cs->msh->x[ip];
		    double x0 = 2.0*cs->msh->x_max[0]/8.0; 		
		    double y  = cs->msh->y[ip]; 		
		    double y0 = 1.2*cs->msh->x_max[1]/2.0; 		
		    double z  = cs->msh->z[ip]; 		
		    double z0 = cs->msh->x_max[2]/2.0; 		

		    double r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0);

                    cs->rho0[ip] = 1.0;//0.4 + 3.0*(x/20.0);
                    cs->p0[ip]   = 1.0/cs->ig->gamma;
                    cs->U0[ip]   = 0.3;
                    cs->V0[ip]   = 0.0;
                    cs->W0[ip]   = 0.0;
                }
            }
        }
    }else{

	string filename = opt->filename;

	MPI_File fh;
	MPI_Offset disp, filesize;

	MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

	//Need to displace by the first 3 doubles and then three double fields
	disp = 5*sizeof(double)+ sizeof(double)*3.0*opt->Nx*opt->Ny*opt->Nz;

	IF_RANK0{
	    cout << " " << endl;
	}

	double *rho_in,
	       *rhoU_in,
	       *rhoV_in,
	       *rhoW_in,
	       *rhoE_in;

	cs->c2d->allocY(rho_in);
	cs->c2d->allocY(rhoU_in);
	cs->c2d->allocY(rhoV_in);
	cs->c2d->allocY(rhoW_in);
	cs->c2d->allocY(rhoE_in);

	cs->c2d->readVar(fh, disp, 1, rho_in);
	cs->c2d->readVar(fh, disp, 1, rhoU_in);
	cs->c2d->readVar(fh, disp, 1, rhoV_in);
	cs->c2d->readVar(fh, disp, 1, rhoW_in);
	cs->c2d->readVar(fh, disp, 1, rhoE_in);

	MPI_File_close(&fh);

        FOR_Z_YPEN{
            FOR_Y_YPEN{
                FOR_X_YPEN{

                    int ip = GETMAJIND_YPEN;
		    int jp = GETIND_YPEN;

	            cs->rho0[ip] = rho_in[jp];
	            cs->U0[ip] = cs->ig->solveU(rho_in[jp], rhoU_in[jp]);
	            cs->V0[ip] = cs->ig->solveU(rho_in[jp], rhoV_in[jp]);
	            cs->W0[ip] = cs->ig->solveU(rho_in[jp], rhoW_in[jp]);
	            cs->p0[ip] = cs->ig->solvep(rho_in[jp], rhoE_in[jp], cs->U0[ip], cs->V0[ip], cs->W0[ip]);

		}
	    }
	}

	cs->c2d->deallocXYZ(rho_in);
	cs->c2d->deallocXYZ(rhoU_in);
	cs->c2d->deallocXYZ(rhoV_in);
	cs->c2d->deallocXYZ(rhoW_in);
	cs->c2d->deallocXYZ(rhoE_in);

    }


    cs->setInitialConditions();



    ///////////////////////////
    //Add images to be output//
    ///////////////////////////


    //This is probably bad programming, but we'll downcast the abstract solver class pointer to the
    //solver pointer so we can access the add image function and the solver member we want to print out
    CurvilinearCSolver *cs_downcast = static_cast<CurvilinearCSolver*>(cs);
    cs_downcast->addImageOutput(new PngWriter(50, 1024, 1024, cs_downcast->p, "P", 2, 0.5, PngWriter::BWR));
    cs_downcast->addImageOutput(new PngWriter(50, 1024, 1024, cs_downcast->V, "V", 2, 0.5, -0.4, 0.4, PngWriter::BWR));
    cs_downcast->addImageOutput(new PngWriter(50, 1024, 1024, cs_downcast->U, "U", 2, 0.5, 0.0, 0.65, PngWriter::RAINBOW));
    cs_downcast->addImageOutput(new PngWriter(50, 1024, 1024, cs->varData[3], "VORTMAG", 2, 0.5, PngWriter::RAINBOW));
    cs_downcast->addImageOutput(new PngWriter(50, 1024, 1024, cs->varData[4], "DIL", 2, 0.5, -0.7,0.1, PngWriter::GREYSCALE));


    ////////////////////////////////////////
    //Execute the solver timestepping loop//
    ////////////////////////////////////////
    rk->executeSolverLoop();  

    //Now lets kill MPI
    MPI_Finalize();

    return 0;
}

void CurvilinearCSolver::initialHook(){

    //Lets do dilatation and vorticity as a test...
    double *vortX, *vortY, *vortZ, *vort_mag, *dil;

    //Add to the variable data list, this is clunky?
    varData.push_back(vortX);
    varData.push_back(vortY);
    varData.push_back(vortZ);
    varData.push_back(vort_mag);
    varData.push_back(dil);

    //Now lets run through the list and allocate the data
    for(vector<double*>::iterator itr = varData.begin(); itr != varData.end(); ++itr){
	c2d->allocY(*itr);
	FOR_XYZ_YPEN{
	    (*itr)[ip] = 0.0;
	}
    }


};


void CurvilinearCSolver::fullStepTemporalHook(){

    //Work with these pointers since its easier....
    double *vortX, *vortY, *vortZ, *vort_mag, *dil;
   
    vortX    = varData[0];
    vortY    = varData[1];
    vortZ    = varData[2];
    vort_mag = varData[3];
    dil      = varData[4];

    if(timeStep%25 == 0){

    //Specifiy our input and output vectors for the computeGradient call
    double *vi[] = {U, V, W};
    vector<double*> vecIn(vi, vi+sizeof(vi)/sizeof(vi[0]));

    double *vo[] = {dU1, dU2, dU3, dV1, dV2, dV3, dW1, dW2, dW3};
    vector<double*> vecOut(vo, vo+sizeof(vo)/sizeof(vo[0]));

    computeGradient(vecIn, vecOut);

    //Get each of the cartesian velocity derivative components
    double *dUdx, *dUdy, *dUdz;
    double *dVdx, *dVdy, *dVdz;
    double *dWdx, *dWdy, *dWdz;

    dUdx = tempY1; dUdy = tempY2; dUdz = tempY3;
    dVdx = tempY4; dVdy = tempY5; dVdz = tempY6;
    dWdx = tempY7; dWdy = tempY8; dWdz = tempY9;

    FOR_XYZ_YPEN{

	dUdx[ip] = J11[ip]*dU1[ip] + J21[ip]*dU2[ip] + J31[ip]*dU3[ip];
	dUdy[ip] = J12[ip]*dU1[ip] + J22[ip]*dU2[ip] + J32[ip]*dU3[ip];
	dUdz[ip] = J13[ip]*dU1[ip] + J23[ip]*dU2[ip] + J33[ip]*dU3[ip];

	dVdx[ip] = J11[ip]*dV1[ip] + J21[ip]*dV2[ip] + J31[ip]*dV3[ip];
	dVdy[ip] = J12[ip]*dV1[ip] + J22[ip]*dV2[ip] + J32[ip]*dV3[ip];
	dVdz[ip] = J13[ip]*dV1[ip] + J23[ip]*dV2[ip] + J33[ip]*dV3[ip];

	dWdx[ip] = J11[ip]*dW1[ip] + J21[ip]*dW2[ip] + J31[ip]*dW3[ip];
	dWdy[ip] = J12[ip]*dW1[ip] + J22[ip]*dW2[ip] + J32[ip]*dW3[ip];
	dWdz[ip] = J13[ip]*dW1[ip] + J23[ip]*dW2[ip] + J33[ip]*dW3[ip];

	dil[ip] = dUdx[ip] + dVdy[ip] + dWdz[ip];

	vortX[ip] = dWdy[ip] - dVdz[ip]; 
	vortY[ip] = dUdz[ip] - dWdx[ip]; 
	vortZ[ip] = dVdx[ip] - dUdy[ip]; 

	vort_mag[ip] = sqrt(vortX[ip]*vortX[ip] + vortY[ip]*vortY[ip] + vortZ[ip]*vortZ[ip]);

   }

    int Nx = pySize[0];
    int Ny = pySize[1];
    int Nz = pySize[2];   
 
    getRange(vortX, "VORTX", Nx, Ny, Nz, mpiRank);
    getRange(vortY, "VORTY", Nx, Ny, Nz, mpiRank);
    getRange(vortZ, "VORTZ", Nx, Ny, Nz, mpiRank);
    getRange(vort_mag, "VORT_MAG", Nx, Ny, Nz, mpiRank);
    getRange(dil, "DIL", Nx, Ny, Nz, mpiRank);

    }

};


void CurvilinearCSolver::subStepTemporalHook(){};
void CurvilinearCSolver::preStepBoundaryHook(){};
void CurvilinearCSolver::postStepBoundaryHook(){};

double CurvilinearCSolver::contRHSSource(int ip){return 0.0;};
double CurvilinearCSolver::xmomRHSSource(int ip){return 0.0;};
double CurvilinearCSolver::ymomRHSSource(int ip){return 0.0;};
double CurvilinearCSolver::zmomRHSSource(int ip){return 0.0;};
double CurvilinearCSolver::engyRHSSource(int ip){return 0.0;};









