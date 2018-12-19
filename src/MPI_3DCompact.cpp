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
 

    ////////////////////////////////////
    //Time Stepping info intialization//
    ////////////////////////////////////
    TimeStepping::TimeSteppingType timeSteppingType = TimeStepping::CONST_CFL;
    double CFL       = 0.5;
    int maxTimeStep  = 25000;
    double maxTime   = 3000.0;
    int filterStep   = 1;
    int checkStep    = 1;
    int dumpStep     = 2000;
    int imageStep    = 50;
    TimeStepping *ts = new TimeStepping(timeSteppingType, 
					CFL, 
					maxTimeStep, 
					maxTime, 
					filterStep, 
					checkStep, 
					dumpStep,
					imageStep, 
					mpiRank);


    ///////////////////////////
    //Boundary Condition Info//
    ///////////////////////////

    BC::BCType bcXType = BC::PERIODIC_SOLVE;
    BC::BCType bcYType = BC::DIRICHLET_SOLVE;
    BC::BCType bcZType = BC::PERIODIC_SOLVE;

    BC::BCKind bcX0 = BC::INTERNALLY_PERIODIC;
    BC::BCKind bcX1 = BC::INTERNALLY_PERIODIC;
    BC::BCKind bcY0 = BC::ADIABATIC_WALL;
    BC::BCKind bcY1 = BC::CYL_CURVILINEARSPONGE;
    BC::BCKind bcZ0 = BC::PERIODIC;
    BC::BCKind bcZ1 = BC::PERIODIC;

    double periodicDisp[3][3];
    //x-direction periodic displacement
    periodicDisp[0][0] = 0.0; // in x
    periodicDisp[0][1] = 0.0; // in y
    periodicDisp[0][2] = 0.0; // in z

    //y-direction periodic displacement
    periodicDisp[1][0] = 0.0; // in x
    periodicDisp[1][1] = 0.0; // in y
    periodicDisp[1][2] = 0.0; // in z

    //z-direction periodic displacement
    periodicDisp[2][0] = 0.0; // in x
    periodicDisp[2][1] = 0.0; // in y
    periodicDisp[2][2] = M_PI; // in z


    bool periodicBC[3];
    BC *bc = new BC(bcXType, bcX0, bcX1,
                    bcYType, bcY0, bcY1,
                    bcZType, bcZ0, bcZ1,
		    periodicBC, mpiRank, periodicDisp);

    /////////////////////////
    //Initialize the Domain//
    /////////////////////////
    int    Nx = 165,
           Ny = 110,
           Nz = 48;

    //For curvilinear coordinates these should all correspond to the max xi, eta, and zeta values
    double Lx = 1.0,
           Ly = 1.0,
           Lz = 1.0;
    Domain *d = new Domain(bc, Nx, Ny, Nz, Lx, Ly, Lz, mpiRank);


    /////////////////////////////
    //Initializing Pencil Decomp//
    //////////////////////////////
 
    int pRow = 0, 
	pCol = 0;
    IF_RANK0 cout << endl << " > Initializing the pencil decomposition... " << endl;

    C2Decomp *c2d;
    c2d = new C2Decomp(Nx, Ny, Nz, pRow, pCol, periodicBC);
  
    IF_RANK0 cout << " > Successfully initialized! " << endl;
    IF_RANK0 cout << " > Handing some decomp info back to the domain object... " << endl;
    d->setPencilDecompInfo(c2d);

    //This is info needed to use the macros...
    int pxSize[3], pySize[3], pzSize[3];
    int pxStart[3], pyStart[3], pzStart[3];
    int pxEnd[3], pyEnd[3], pzEnd[3];
    d->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);

    /////////////////////////
    //Initialize the Solver//
    /////////////////////////
    double alphaF  = 0.45;
    double mu_ref  = 0.001;
    bool useTiming = false;
    AbstractCSolver *cs= new CurvilinearCSolver(c2d, d, bc, ts, alphaF, mu_ref, useTiming);

    //Attach the mesh object to the solver...
    cs->msh = new AlgebraicSingleBlockMesh(c2d, cs, d, mpiRank);

    ///////////////////////////////////////////
    //Initialize Execution Loop and RK Method//
    ///////////////////////////////////////////
    AbstractRK *rk = new TVDRK3(cs);

    ///////////////////////////////
    //Set flow initial conditions//
    ///////////////////////////////

    bool fromRestart = false;

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

	string filename = "SolutionDump.2000";
	int timestep_start = 2000;

	cs->timeStep = timestep_start;

	MPI_File fh;
	MPI_Offset disp, filesize;

	MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	disp = 0;

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

	//Segfaults for some reason if I dealloc these?
	//cs->c2d->deallocXYZ(rhoW_in);
	//cs->c2d->deallocXYZ(rhoE_in);

    }


    cs->setInitialConditions();



    ///////////////////////////
    //Add images to be output//
    ///////////////////////////


    //This is probably bad programming, but we'll downcast the abstract solver class pointer to the
    //solver pointer so we can access the add image function and the solver member we want to print out
    CurvilinearCSolver *cs_downcast = static_cast<CurvilinearCSolver*>(cs);
    cs_downcast->addImageOutput(new PngWriter(25, 1028, 1028, cs_downcast->p, "P", 2, 0.5, PngWriter::BWR));
    cs_downcast->addImageOutput(new PngWriter(25, 1028, 1028, cs_downcast->V, "V", 2, 0.5, -0.4, 0.4, PngWriter::BWR));
    cs_downcast->addImageOutput(new PngWriter(25, 1028, 1028, cs_downcast->U, "U", 2, 0.5, 0.0, 0.65, PngWriter::RAINBOW));
    cs_downcast->addImageOutput(new PngWriter(25, 1028, 1028, cs->varData[3], "VORTMAG", 2, 0.5, PngWriter::RAINBOW));
    cs_downcast->addImageOutput(new PngWriter(25, 1028, 1028, cs->varData[4], "DIL", 2, 0.5, -2.0,2.0, PngWriter::BWR));


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
    
    /////////////////////////////
    //Xi2-Direction Derivatives//
    /////////////////////////////

    //First we'll do all of the ~Y-Direction derivatives to calc tau
    derivXi2->calc1stDerivField(U, dU2);
    derivXi2->calc1stDerivField(V, dV2);
    derivXi2->calc1stDerivField(W, dW2);

    /////////////////////////////
    //Xi1-Direction Derivatives//
    /////////////////////////////

    double *dU1_xp, *dV1_xp, *dW1_xp;

    //Point to the needed X memory
    U_xp = tempX1; dU1_xp = tempX5;
    V_xp = tempX2; dV1_xp = tempX6;
    W_xp = tempX3; dW1_xp = tempX7;

    c2d->transposeY2X_MajorIndex(U, U_xp);
    c2d->transposeY2X_MajorIndex(V, V_xp);
    c2d->transposeY2X_MajorIndex(W, W_xp);

    derivXi1->calc1stDerivField(U_xp, dU1_xp);
    derivXi1->calc1stDerivField(V_xp, dV1_xp);
    derivXi1->calc1stDerivField(W_xp, dW1_xp);

    c2d->transposeX2Y_MajorIndex(dU1_xp, dU1);
    c2d->transposeX2Y_MajorIndex(dV1_xp, dV1);
    c2d->transposeX2Y_MajorIndex(dW1_xp, dW1);

    /////////////////////////////
    //Xi3-Direction Derivatives//
    /////////////////////////////

    double *dU3_zp, *dV3_zp, *dW3_zp;

    //Point to the needed Z memory
    U_zp = tempZ1; dU3_zp = tempZ5;
    V_zp = tempZ2; dV3_zp = tempZ6;
    W_zp = tempZ3; dW3_zp = tempZ7;

    c2d->transposeY2Z_MajorIndex(U, U_zp);
    c2d->transposeY2Z_MajorIndex(V, V_zp);
    c2d->transposeY2Z_MajorIndex(W, W_zp);

    derivXi3->calc1stDerivField(U_zp, dU3_zp);
    derivXi3->calc1stDerivField(V_zp, dV3_zp);
    derivXi3->calc1stDerivField(W_zp, dW3_zp);

    c2d->transposeZ2Y_MajorIndex(dU3_zp, dU3);
    c2d->transposeZ2Y_MajorIndex(dV3_zp, dV3);
    c2d->transposeZ2Y_MajorIndex(dW3_zp, dW3);


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









