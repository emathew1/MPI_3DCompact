#include "SpongeBC.hpp"

SpongeBC::SpongeBC(AbstractSingleBlockMesh *msh, Domain *domain, IdealGas *idealGas, BC *bc, C2Decomp *c2d, int mpiRank){

	    
	    IF_RANK0 std::cout << endl;
	    IF_RANK0 std::cout << " > Sponge BC found, initializing Sponge average fields and strength fields..." << std::endl;
	

	    this->msh = msh;
	    this->domain = domain;
	    this->idealGas = idealGas;
	    this->bc = bc;

	    domain->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);

	    

	    this->Nx = domain->gNx;
	    this->Ny = domain->gNy;
	    this->Nz = domain->gNz;
	    N = Nx*Ny*Nz;

	    this->mpiRank = mpiRank;

	    c2d->allocY(sigma);
	    c2d->allocY(spongeRhoAvg);
	    c2d->allocY(spongeRhoUAvg);
	    c2d->allocY(spongeRhoVAvg);
	    c2d->allocY(spongeRhoWAvg);
	    c2d->allocY(spongeRhoEAvg);

	    avgT = 1.0;
	    epsP = 0.005;
	    spongeP = 1.0/idealGas->gamma;
	    spongeStrength = 2;

	    //Need to initialize the sponge sigma to zero
	    FOR_XYZ_YPEN sigma[ip] = 0.0;

 	    //Need to have some way to trigger which one of these is inialized...

	    //If rectangular sponge BC
	    initRectSpongeBC();
 
	    //If cylindrical sponge BC
	    initCylSpongeBC();	   

	    //If spherical sponge BC


            getRange(sigma, "SPONGE SIGMA", pySize[0], pySize[1], pySize[2], mpiRank);
	    IF_RANK0 std::cout << " > Done initializing sponge!" << std::endl;

	}


void SpongeBC::initRectSpongeBC(){


	    //Default the maximum ends of the sponge to the domain max, can and may be changed for curvilinear domains
	    double spongeXMin = msh->x_min[0];
	    double spongeXMax = msh->x_max[0];

	    double spongeYMin = msh->x_min[1];
	    double spongeYMax = msh->x_max[1];

	    double spongeZMin = msh->x_min[2];
	    double spongeZMax = msh->x_max[2];

	    //Default to an 1/8th of the domain size in that direction 
	    spongeLX = 0.125*(spongeXMax-spongeXMin);
	    spongeLY = 0.125*(spongeYMax-spongeYMin);
	    spongeLZ = 0.125*(spongeZMax-spongeZMin);

	    IF_RANK0 cout << " > spongeLX = " << spongeLX << endl;	    
	    IF_RANK0 cout << " > spongeLY = " << spongeLY << endl;	    
	    IF_RANK0 cout << " > spongeLZ = " << spongeLZ << endl;	    
	
	    //Use this data to initialize the sponge zones / sponge sigma strength...
	    if(bc->bcX0 == Options::SPONGE){
		FOR_X_YPEN{
		    FOR_Y_YPEN{
			FOR_Z_YPEN{
		    	    int ip = GETMAJIND_YPEN;
			    double dx = msh->x[ip]-spongeXMin;
		            if(dx < spongeLX){
		        	double spongeX = (spongeLX - dx)/spongeLX;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ip]);
			    }
		        }
		    }
	  	}
	    }

	    if(bc->bcX1 == Options::SPONGE){
	        FOR_X_YPEN{
		    FOR_Y_YPEN{
		        FOR_Z_YPEN{
			    int ip = GETMAJIND_YPEN;
			    double dx = spongeXMax - msh->x[ip];
		    	    if(dx < spongeLX){
		        	double spongeX = (msh->x[ip] - (spongeXMax - spongeLX))/spongeLX;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ip]);
			    }
		        }
		    }
	  	}
	    }     

	    if(bc->bcY0 == Options::SPONGE){
	        FOR_X_YPEN{
	            FOR_Y_YPEN{
			FOR_Z_YPEN{
		            int ip = GETMAJIND_YPEN;
			    double dy = msh->y[ip]-spongeYMin;	
		            if(dy < spongeLY){
		                double spongeY = (spongeLY - dy)/spongeLY;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ip]);
			    }
		        }
		    }
		}
	    }
	
	    if(bc->bcY1 == Options::SPONGE){
		FOR_X_YPEN{
		    FOR_Y_YPEN{
		        FOR_Z_YPEN{
		            int ip = GETMAJIND_YPEN;
			    double dy = spongeYMax-msh->y[ip];
		            if(dy < spongeLY){
		                double spongeY = (msh->y[ip] - (spongeYMax - spongeLY))/spongeLY;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ip]);
			    }
		        }
		    }
		}
	    }    

	    if(bc->bcZ0 == Options::SPONGE){
		FOR_X_YPEN{
		    FOR_Y_YPEN{
		        FOR_Z_YPEN{
		            int ip = GETMAJIND_YPEN;
			    double dz = msh->z[ip] - spongeZMin;
		            if(dz < spongeLZ){
		                double spongeZ = (spongeLZ - dz)/spongeLZ;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeZ, 2.0) + 0.845*pow(spongeZ, 8.0)), sigma[ip]);
			    }
		        }
		    }
		}
  	    }
	
	    if(bc->bcZ1 == Options::SPONGE){
	        FOR_X_YPEN{
	            FOR_Y_YPEN{
		        FOR_Z_YPEN{
		    	    int ip = GETMAJIND_YPEN;
			    double dz = spongeZMax-msh->z[ip];
		    	    if(dz < spongeLZ){
		        	double spongeZ = (msh->z[ip] - (spongeZMax-spongeLZ))/spongeLZ;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeZ, 2.0) + 0.845*pow(spongeZ, 8.0)), sigma[ip]);
			    }
		        }
		    }
		}
	    }    


}

void SpongeBC::initCylSpongeBC(){


	    //Default the maximum ends of the sponge to the domain max, can and may be changed for curvilinear domains
	    double spongeXMin = msh->x_min[0];
	    double spongeXMax = msh->x_max[0];
	    double spongeXAvg = 0.5*(spongeXMin + spongeXMax);
	    double spongeXRmax = spongeXMax;
	    double spongeXRmin = 5.0;

	    double spongeYMin = msh->x_min[1];
	    double spongeYMax = msh->x_max[1];
	    double spongeYAvg = 0.5*(spongeYMin + spongeYMax);
	    double spongeYRmax = spongeYMax;
	    double spongeYRmin = 5.0;

	    double spongeZMin = msh->x_min[2];
	    double spongeZMax = msh->x_max[2];
	    double spongeZAvg = 0.5*(spongeZMin + spongeZMax);
	    double spongeZRmax = spongeZMax;
	    double spongeZRmin = 5.0;
		
	    if(bc->bcY1 == Options::SPONGE){
		FOR_X_YPEN{
		    FOR_Y_YPEN{
			FOR_Z_YPEN{
			    int ip = GETMAJIND_YPEN;
			    double r = sqrt(pow(msh->x[ip]-spongeXAvg,2.0) + pow(msh->y[ip]-spongeYAvg,2.0));
			    if(r > spongeYRmin){
				double spongeR = (r-spongeYRmin)/(spongeYRmax-spongeYRmin);
				sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeR,2.0) + 0.845*pow(spongeR, 8.0)), sigma[ip]);
			    }
			}
		    }
		}
	    }

	    //Default to an 1/8th of the domain size in that direction 
/*	
	    //Use this data to initialize the sponge zones / sponge sigma strength...
	    if(bc->bcX0 == Options::RECT_CURVILINEARSPONGE){
		FOR_X_YPEN{
		    FOR_Y_YPEN{
			FOR_Z_YPEN{
		    	    int ip = GETMAJIND_YPEN;
			    double dx = msh->x[ip]-spongeXMin;
		            if(dx < spongeLX){
		        	double spongeX = (spongeLX - dx)/spongeLX;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ip]);
			    }
		        }
		    }
	  	}
	    }

	    if(bc->bcX1 == Options::RECT_CURVILINEARSPONGE){
	        FOR_X_YPEN{
		    FOR_Y_YPEN{
		        FOR_Z_YPEN{
			    int ip = GETMAJIND_YPEN;
			    double dx = spongeXMax - msh->x[ip];
		    	    if(dx < spongeLX){
		        	double spongeX = (msh->x[ip] - (spongeXMax - spongeLX))/spongeLX;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ip]);
			    }
		        }
		    }
	  	}
	    }     

	    if(bc->bcY0 == Options::RECT_CURVILINEARSPONGE){
	        FOR_X_YPEN{
	            FOR_Y_YPEN{
			FOR_Z_YPEN{
		            int ip = GETMAJIND_YPEN;
			    double dy = msh->y[ip]-spongeYMin;	
		            if(dy < spongeLY){
		                double spongeY = (spongeLY - dy)/spongeLY;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ip]);
			    }
		        }
		    }
		}
	    }
	
	    if(bc->bcY1 == Options::RECT_CURVILINEARSPONGE){
		FOR_X_YPEN{
		    FOR_Y_YPEN{
		        FOR_Z_YPEN{
		            int ip = GETMAJIND_YPEN;
			    double dy = spongeYMax-msh->y[ip];
		            if(dy < spongeLY){
		                double spongeY = (msh->y[ip] - (spongeYMax - spongeLY))/spongeLY;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ip]);
			    }
		        }
		    }
		}
	    }    

	    if(bc->bcZ0 == Options::RECT_CURVILINEARSPONGE){
		FOR_X_YPEN{
		    FOR_Y_YPEN{
		        FOR_Z_YPEN{
		            int ip = GETMAJIND_YPEN;
			    double dz = msh->z[ip] - spongeZMin;
		            if(dz < spongeLZ){
		                double spongeZ = (spongeLZ - dz)/spongeLZ;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeZ, 2.0) + 0.845*pow(spongeZ, 8.0)), sigma[ip]);
			    }
		        }
		    }
		}
  	    }
	
	    if(bc->bcZ1 == Options::RECT_CURVILINEARSPONGE){
	        FOR_X_YPEN{
	            FOR_Y_YPEN{
		        FOR_Z_YPEN{
		    	    int ip = GETMAJIND_YPEN;
			    double dz = spongeZMax-msh->z[ip];
		    	    if(dz < spongeLZ){
		        	double spongeZ = (msh->z[ip] - (spongeZMax-spongeLZ))/spongeLZ;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeZ, 2.0) + 0.845*pow(spongeZ, 8.0)), sigma[ip]);
			    }
		        }
		    }
		}
	    }    
*/

	}




