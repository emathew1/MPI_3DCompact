#ifndef STATS_HPP
#define STATS_HPP

#include "Macros.hpp"
#include "Domain.hpp"
#include "Options.hpp"
#include "CurvilinearCSolver.hpp"

class Stats{

    public:

	CurvilinearCSolver *cs;
	Options *opt;

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

	Stats(CurvilinearCSolver *cs, Options *opt){

	    this->cs  = cs;
	    this->opt = opt;

	    //Get local copies of the domain sizes
            cs->dom->getPencilDecompInfo(pxSize, pySize, pzSize, pxStart, pyStart, pzStart, pxEnd, pyEnd, pzEnd);



	    //Get the links to the correct source info
	    U = cs->U; 
	    V = cs->V; 
	    W = cs->W; 
 
	    rho = cs->rho;
	    p   = cs->p;
	    T   = cs->T;

	    //get the flags for the statistics
	    velocityStats = opt->velocityStats;
	    thermoStats   = opt->thermoStats;

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

		}else{
		
		    thermoStatsWeight = 0.0;
		    //initialize the fields
		    RHOAVG[ip] = 0.0;
		    RHORMS[ip] = 0.0;
		    PAVG[ip] = 0.0;
		    PRMS[ip] = 0.0;
		    TAVG[ip] = 0.0;
		    TRMS[ip] = 0.0;
		}

	    }
	    
	}


	void dumpStatsFields();
	void updateStatsFields();

};

#endif
