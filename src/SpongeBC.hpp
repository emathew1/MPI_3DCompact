#ifndef _SPONGEBCH_
#define _SPONGEBCH_

#include <iostream>
#include "Macros.hpp"
#include "Domain.hpp"
#include "IdealGas.hpp"
#include "BC.hpp"

class SpongeBC{

    public:

	Domain *domain;
	IdealGas *idealGas;
	BC *bc;

	int Nx, Ny, Nz, N;

        int pxSize[3], pySize[3], pzSize[3]; 
        int pxStart[3], pyStart[3], pzStart[3];
        int pxEnd[3], pyEnd[3], pzEnd[3];

	double avgT;
	double epsP;
	double spongeP;
	double spongeStrength;
	double spongeLX;
	double spongeLY;
	double spongeLZ;

	double *sigma;
	double *spongeRhoAvg;
	double *spongeRhoUAvg;
	double *spongeRhoVAvg;
	double *spongeRhoWAvg;
	double *spongeRhoEAvg;
    
	SpongeBC(Domain *domain, IdealGas *idealGas, BC *bc, C2Decomp *c2d, int mpiRank){

	    
	    IF_RANK0 std::cout << endl;
	    IF_RANK0 std::cout << " > Sponge BC found, initializing Sponge average fields and strength fields..." << std::endl;
	
	    this->domain = domain;
	    this->idealGas = idealGas;
	    this->bc = bc;

	    domain->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);


	    this->Nx = domain->gNx;
	    this->Ny = domain->gNy;
	    this->Nz = domain->gNz;
	    N = Nx*Ny*Nz;

	    c2d->allocX(sigma);
	    c2d->allocX(spongeRhoAvg);
	    c2d->allocX(spongeRhoUAvg);
	    c2d->allocX(spongeRhoVAvg);
	    c2d->allocX(spongeRhoWAvg);
	    c2d->allocX(spongeRhoEAvg);

	    avgT = 15.0;
	    epsP = 0.005;
	    spongeP = 1.0/idealGas->gamma;
	    spongeStrength = 6.0;
	    spongeLX = 0.1*domain->gLx;
	    spongeLY = 0.1*domain->gLy;
	    spongeLZ = 0.1*domain->gLz;

	    //Need to initialize the sponge sigma to zero
	    FOR_XYZ_XPEN sigma[ip] = 0.0;
	
	    //Use this data to initialize the sponge zones / sponge sigma strength...
	    if(bc->bcX0 == BC::SPONGE){
		FOR_X_XPEN{
		    int ip = GETGLOBALXIND_XPEN;
		    if(domain->x[ip] < spongeLX){
		        double spongeX = (spongeLX - domain->x[ip])/spongeLX;
		        FOR_Y_XPEN{
			    FOR_Z_XPEN{
			        int ii = GETMAJIND_XPEN;
			        sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ii]);
			    }
		        }
		    }
		}
	    }

	    if(bc->bcX1 == BC::SPONGE){
		FOR_X_XPEN{
		    int ip = GETGLOBALXIND_XPEN;
		    if(domain->x[ip] > domain->gLx - spongeLX){
		        double spongeX = (domain->x[ip] - (domain->gLx - spongeLX))/spongeLX;
		        FOR_Y_XPEN{
			    FOR_Z_XPEN{
			        int ii = GETMAJIND_XPEN;
			        sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ii]);
			    }
		        }
		    }
		}
	    }    

	    if(bc->bcY0 == BC::SPONGE){
		FOR_Y_XPEN{
		    int jp = GETGLOBALYIND_XPEN;
		    if(domain->y[jp] < spongeLY){
		        double spongeY = (spongeLY - domain->y[jp])/spongeLY;
		        FOR_X_XPEN{
			    FOR_Z_XPEN{
			        int ii = GETMAJIND_XPEN;
			        sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ii]);
			    }
		        }
		    }
		}
	    }
	
	    if(bc->bcY1 == BC::SPONGE){
		FOR_Y_XPEN{
		    int jp = GETGLOBALYIND_XPEN;
		    if(domain->y[jp] > domain->gLy - spongeLY){
		        double spongeY = (domain->y[jp] - (domain->gLy - spongeLY))/spongeLY;
		        FOR_X_XPEN{
			    FOR_Z_XPEN{
			        int ii = GETMAJIND_XPEN;
			        sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ii]);
			    }
		        }
		    }
		}
	    }    

	    if(bc->bcZ0 == BC::SPONGE){
		FOR_Z_XPEN{
		    int kp = GETGLOBALZIND_XPEN;
		    if(domain->z[kp] < spongeLZ){
		        double spongeZ = (spongeLZ - domain->z[kp])/spongeLZ;
		        FOR_X_XPEN{
			    FOR_Y_XPEN{
			        int ii = GETMAJIND_XPEN;
			        sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeZ, 2.0) + 0.845*pow(spongeZ, 8.0)), sigma[ii]);
			    }
		        }
		    }
		}
	    }
	
	    if(bc->bcZ1 == BC::SPONGE){
		FOR_Z_XPEN{
		    int kp = GETGLOBALZIND_XPEN;
		    if(domain->z[kp] > domain->gLz - spongeLZ){
		        double spongeZ = (domain->z[kp] - (domain->gLz - spongeLZ))/spongeLZ;
		        FOR_X_XPEN{
			    FOR_Y_XPEN{
			        int ii = GETMAJIND_XPEN;
			        sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeZ, 2.0) + 0.845*pow(spongeZ, 8.0)), sigma[ii]);
			    }
		        }
		    }
		}
	    }    



	    IF_RANK0 std::cout << " > Done initializing sponge!" << std::endl;

	}

};

#endif
