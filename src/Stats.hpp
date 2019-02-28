#ifndef STATS_HPP
#define STATS_HPP

#include "Macros.hpp"
#include "Domain.hpp"
#include "Options.hpp"
#include "AbstractCSolver.hpp"

class Stats{

    public:

	AbstractCSolver *cs;
	Options *opt;

	int mpiRank;

        int pxSize[3], pySize[3], pzSize[3];
        int pxStart[3], pyStart[3], pzStart[3];
        int pxEnd[3], pyEnd[3], pzEnd[3];

	bool velocityStats;
	double velocityStatsWeight;
	string velocityStatsFilename;
	double *U, *V, *W;
	double *UAVG, *VAVG, *WAVG;
	double *URMS, *VRMS, *WRMS;
	double *UVREY, *UWREY, *VWREY;

	bool thermoStats;
	double thermoStatsWeight;
	string thermoStatsFilename;
	double *rho, *p, *T;
	double *RHOAVG, *RHORMS;
	double *PAVG, *PRMS;
	double *TAVG, *TRMS;

	Options::StatsAvgType statsAvgType; 
	int stats_interval; 
	double prev_time;

	Stats(AbstractCSolver *cs, Options *opt){

	    this->cs  = cs;
	    this->opt = opt;

	    //Get local copies of the domain sizes
            cs->dom->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);

	    mpiRank = cs->mpiRank;

	    //Get the links to the correct source info
	    U = cs->U; 
	    V = cs->V; 
	    W = cs->W; 
 
	    rho = cs->rho2;
	    p   = cs->p;
	    T   = cs->T;

	    //get the flags for the statistics
	    velocityStats = opt->velocityStats;
	    thermoStats   = opt->thermoStats;

	    //get the filenames for the statistics
	    velocityStatsFilename = opt->velocityStatsFilename;
	    thermoStatsFilename = opt->thermoStatsFilename;

	    //Initialize the previous time to the current simulation time...
	    prev_time = cs->time;

	    //If from restart file
	    if(opt->velocityStats){

		//allocate the arrays...
		cs->c2d->allocY(UAVG);	
		cs->c2d->allocY(URMS);	
		cs->c2d->allocY(VAVG);	
		cs->c2d->allocY(VRMS);	
		cs->c2d->allocY(WAVG);	
		cs->c2d->allocY(WRMS);	
		cs->c2d->allocY(UVREY);	
		cs->c2d->allocY(UWREY);	
		cs->c2d->allocY(VWREY);	

	        if(opt->velocityStatsFromRestart){
		     //read in the weight and the stat fields from the files
		    IF_RANK0{
			FILE *ptr;
			ptr = fopen(velocityStatsFilename.c_str(), "rb");
			if(ptr == NULL){
			    cout << "ERROR: Couldn't open file " << velocityStatsFilename << endl;
			    MPI_Abort(MPI_COMM_WORLD, -10);

			}else{
			    fread(&velocityStatsWeight, sizeof(double), 1, ptr);
			}
			fclose(ptr);
		    }
		    MPI_Bcast(&velocityStatsWeight, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


		    MPI_File fh;
		    MPI_Offset disp, filesize;
		    MPI_File_open(MPI_COMM_WORLD, velocityStatsFilename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

		    disp = sizeof(double);
	
		    cs->c2d->readVar(fh, disp, 1, cs->tempY1); 
		    cs->c2d->readVar(fh, disp, 1, cs->tempY2); 
		    cs->c2d->readVar(fh, disp, 1, cs->tempY3); 
		    cs->c2d->readVar(fh, disp, 1, cs->tempY4); 
		    cs->c2d->readVar(fh, disp, 1, cs->tempY5); 
		    cs->c2d->readVar(fh, disp, 1, cs->tempY6); 
		    cs->c2d->readVar(fh, disp, 1, cs->tempY7); 
		    cs->c2d->readVar(fh, disp, 1, cs->tempY8); 
		    cs->c2d->readVar(fh, disp, 1, cs->tempY9); 
		
		    MPI_File_close(&fh);

		    FOR_Z_YPEN{
			FOR_Y_YPEN{
			    FOR_X_YPEN{
				int ip = GETMAJIND_YPEN;
				int jp = GETIND_YPEN;
				
				UAVG[ip] = cs->tempY1[jp];
				VAVG[ip] = cs->tempY2[jp];
				WAVG[ip] = cs->tempY3[jp];
				URMS[ip] = cs->tempY4[jp];
				VRMS[ip] = cs->tempY5[jp];
				WRMS[ip] = cs->tempY6[jp];
				UVREY[ip] = cs->tempY7[jp];
				UWREY[ip] = cs->tempY8[jp];
				VWREY[ip] = cs->tempY9[jp];

			    }
			}
		    }

	        }else{
		
		    velocityStatsWeight = 0.0;
		    //initialize the fields
		    FOR_XYZ_YPEN{
			UAVG[ip] = 0.0;
			URMS[ip] = 0.0;
			VAVG[ip] = 0.0;
			VRMS[ip] = 0.0;
			WAVG[ip] = 0.0;
			WRMS[ip] = 0.0;
		        UVREY[ip] = 0.0;	
		        UWREY[ip] = 0.0;	
		        VWREY[ip] = 0.0;	
		    }		    
		}
	    }	

	    IF_RANK0{
		if(opt->velocityStats){
		    cout << " > Initialized Velocity Stats, Weight = " << velocityStatsWeight << endl;
		    if(opt->velocityStatsFromRestart){
		        cout << " > Initialized from file: " << velocityStatsFilename << endl;
		    }
		}
	    }

	    if(opt->thermoStats){

		//allocate the arrays
		cs->c2d->allocY(RHOAVG);	
		cs->c2d->allocY(RHORMS);	
		cs->c2d->allocY(PAVG);	
		cs->c2d->allocY(PRMS);	
		cs->c2d->allocY(TAVG);	
		cs->c2d->allocY(TRMS);	

		//allocate the arrays...
		if(opt->thermoStatsFromRestart){
		    //read in the weight and the stat fields from the files

		    IF_RANK0{
			FILE *ptr;
			ptr = fopen(thermoStatsFilename.c_str(), "rb");
			if(ptr == NULL){
			    cout << "ERROR: Couldn't open file " << thermoStatsFilename << endl;
			    MPI_Abort(MPI_COMM_WORLD, -10);

			}else{
			    fread(&thermoStatsWeight, sizeof(double), 1, ptr);
			}
			fclose(ptr);
		    }
		    MPI_Bcast(&thermoStatsWeight, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


		    MPI_File fh;
		    MPI_Offset disp, filesize;
		    MPI_File_open(MPI_COMM_WORLD, thermoStatsFilename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

		    disp = sizeof(double);
	
		    cs->c2d->readVar(fh, disp, 1, cs->tempY1); 
		    cs->c2d->readVar(fh, disp, 1, cs->tempY2); 
		    cs->c2d->readVar(fh, disp, 1, cs->tempY3); 
		    cs->c2d->readVar(fh, disp, 1, cs->tempY4); 
		    cs->c2d->readVar(fh, disp, 1, cs->tempY5); 
		    cs->c2d->readVar(fh, disp, 1, cs->tempY6); 
		
		    MPI_File_close(&fh);

		    FOR_Z_YPEN{
			FOR_Y_YPEN{
			    FOR_X_YPEN{
				int ip = GETMAJIND_YPEN;
				int jp = GETIND_YPEN;
				
				RHOAVG[ip] = cs->tempY1[jp];
				RHORMS[ip] = cs->tempY2[jp];
				PAVG[ip] = cs->tempY3[jp];
				PRMS[ip] = cs->tempY4[jp];
				TAVG[ip] = cs->tempY5[jp];
				TRMS[ip] = cs->tempY6[jp];

			    }
			}
		    }

		}else{
		
		    thermoStatsWeight = 0.0;
		    //initialize the fields
		    FOR_XYZ_YPEN{
		        RHOAVG[ip] = 0.0;
		        RHORMS[ip] = 0.0;
		        PAVG[ip] = 0.0;
		        PRMS[ip] = 0.0;
		        TAVG[ip] = 0.0;
		        TRMS[ip] = 0.0;
		    }
		}

	    }
	

	    IF_RANK0{
		if(opt->thermoStats){
		    cout << " > Initialized Thermo Stats, Weight = " << thermoStatsWeight << endl;
		    if(opt->thermoStatsFromRestart){
		        cout << " > Initialized from file: " << thermoStatsFilename << endl;
		    }
		}
	    }

	}


	void dumpStatsFields();
	void updateStatsFields();
	void doAveraging(double *phi);
};

#endif
