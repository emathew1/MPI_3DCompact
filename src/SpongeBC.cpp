#include "SpongeBC.hpp"

SpongeBC::SpongeBC(Domain *domain, IdealGas *idealGas, BC *bc, C2Decomp *c2d, int baseDirection, int mpiRank){


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

            if(baseDirection == 0){
                c2d->allocX(sigma);
                c2d->allocX(spongeRhoAvg);
                c2d->allocX(spongeRhoUAvg);
                c2d->allocX(spongeRhoVAvg);
                c2d->allocX(spongeRhoWAvg);
                c2d->allocX(spongeRhoEAvg);
            }else if(baseDirection == 1){
                c2d->allocY(sigma);
                c2d->allocY(spongeRhoAvg);
                c2d->allocY(spongeRhoUAvg);
                c2d->allocY(spongeRhoVAvg);
                c2d->allocY(spongeRhoWAvg);
                c2d->allocY(spongeRhoEAvg);
            }else if(baseDirection == 2){
                c2d->allocZ(sigma);
                c2d->allocZ(spongeRhoAvg);
                c2d->allocZ(spongeRhoUAvg);
                c2d->allocZ(spongeRhoVAvg);
                c2d->allocZ(spongeRhoWAvg);
                c2d->allocZ(spongeRhoEAvg);
            }else{
                cout << "Unknown baseDirection in spongeBC constructor!" << endl;
            }

            avgT = 10.0;
            epsP = 0.005;
            spongeP = 1.0/idealGas->gamma;
            spongeStrength = 12.0;
            spongeLX = 0.125*domain->gLx;
            spongeLY = 0.125*domain->gLy;
            spongeLZ = 0.125*domain->gLz;

            if(baseDirection == 0){

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

            }else if(baseDirection == 1){

                //Need to initialize the sponge sigma to zero
                FOR_XYZ_YPEN sigma[ip] = 0.0;

                //Use this data to initialize the sponge zones / sponge sigma strength...
                if(bc->bcX0 == BC::SPONGE){
                    FOR_X_YPEN{
                        int ip = GETGLOBALXIND_YPEN;
                        if(domain->x[ip] < spongeLX){
                            double spongeX = (spongeLX - domain->x[ip])/spongeLX;
                            FOR_Y_YPEN{
                                FOR_Z_YPEN{
                                    int ii = GETMAJIND_YPEN;
                                    sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ii]);
                                }
                            }
                        }
                    }
                }

                if(bc->bcX1 == BC::SPONGE){
                    FOR_X_YPEN{
                        int ip = GETGLOBALXIND_YPEN;
                        if(domain->x[ip] > domain->gLx - spongeLX){
                            double spongeX = (domain->x[ip] - (domain->gLx - spongeLX))/spongeLX;
                            FOR_Y_YPEN{
                                FOR_Z_YPEN{
                                    int ii = GETMAJIND_YPEN;
                                    sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ii]);
                                }
                            }
                        }
                    }
                }

                if(bc->bcY0 == BC::SPONGE){
                    FOR_Y_YPEN{
                        int jp = GETGLOBALYIND_YPEN;
                        if(domain->y[jp] < spongeLY){
                            double spongeY = (spongeLY - domain->y[jp])/spongeLY;
                            FOR_X_YPEN{
                                FOR_Z_YPEN{
                                    int ii = GETMAJIND_YPEN;
                                    sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ii]);
                                }
                            }
                        }
                    }
                }

                if(bc->bcY1 == BC::SPONGE){
                    FOR_Y_YPEN{
                        int jp = GETGLOBALYIND_YPEN;
                        if(domain->y[jp] > domain->gLy - spongeLY){
                            double spongeY = (domain->y[jp] - (domain->gLy - spongeLY))/spongeLY;
                            FOR_X_YPEN{
                                FOR_Z_YPEN{
                                    int ii = GETMAJIND_YPEN;
                                    sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ii]);
                                }
                            }
                        }
                    }
                }

                if(bc->bcZ0 == BC::SPONGE){
                    FOR_Z_YPEN{
                        int kp = GETGLOBALZIND_YPEN;
                        if(domain->z[kp] < spongeLZ){
                            double spongeZ = (spongeLZ - domain->z[kp])/spongeLZ;
                            FOR_X_YPEN{
                                FOR_Y_YPEN{
                                    int ii = GETMAJIND_YPEN;
                                    sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeZ, 2.0) + 0.845*pow(spongeZ, 8.0)), sigma[ii]);
                                }
                            }
                        }
                    }
                }

                if(bc->bcZ1 == BC::SPONGE){
                    FOR_Z_YPEN{
                        int kp = GETGLOBALZIND_YPEN;
                        if(domain->z[kp] > domain->gLz - spongeLZ){
                            double spongeZ = (domain->z[kp] - (domain->gLz - spongeLZ))/spongeLZ;
                            FOR_X_YPEN{
                                FOR_Y_YPEN{
                                    int ii = GETMAJIND_YPEN;
                                    sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeZ, 2.0) + 0.845*pow(spongeZ, 8.0)), sigma[ii]);
                                }
                            }
                        }
                    }
                }



            }

            IF_RANK0 std::cout << " > Done initializing sponge!" << std::endl;

}




CurvilinearSpongeBC::CurvilinearSpongeBC(AbstractSingleBlockMesh *msh, Domain *domain, IdealGas *idealGas, BC *bc, C2Decomp *c2d, int mpiRank){

	    
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


void CurvilinearSpongeBC::initRectSpongeBC(){


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
	    if(bc->bcX0 == BC::RECT_CURVILINEARSPONGE){
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

	    if(bc->bcX1 == BC::RECT_CURVILINEARSPONGE){
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

	    if(bc->bcY0 == BC::RECT_CURVILINEARSPONGE){
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
	
	    if(bc->bcY1 == BC::RECT_CURVILINEARSPONGE){
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

	    if(bc->bcZ0 == BC::RECT_CURVILINEARSPONGE){
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
	
	    if(bc->bcZ1 == BC::RECT_CURVILINEARSPONGE){
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

void CurvilinearSpongeBC::initCylSpongeBC(){


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
		
	    if(bc->bcY1 == BC::CYL_CURVILINEARSPONGE){
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
	    if(bc->bcX0 == BC::RECT_CURVILINEARSPONGE){
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

	    if(bc->bcX1 == BC::RECT_CURVILINEARSPONGE){
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

	    if(bc->bcY0 == BC::RECT_CURVILINEARSPONGE){
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
	
	    if(bc->bcY1 == BC::RECT_CURVILINEARSPONGE){
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

	    if(bc->bcZ0 == BC::RECT_CURVILINEARSPONGE){
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
	
	    if(bc->bcZ1 == BC::RECT_CURVILINEARSPONGE){
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




