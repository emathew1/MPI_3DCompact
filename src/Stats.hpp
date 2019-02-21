#ifndef STATS_HPP
#define STATS_HPP

#include "Options.hpp"
#include "CurvilinearCSolver.hpp"

class Stats{

    public:

	CurvilinearCSolver *cs;
	Options *opt;


	bool velocityStats;
	double velocityStatsWeight;
	string velocityStatsFilename;
	double *UAVG, *VAVG, *WAVG;
	double *URMS, *VRMS, *WRMS;
	double *UVREY, *UWREY, *VWREY;

	bool thermoStats;
	double thermoStatsWeight;
	string thermoStatsFilename;
	double *RHOAVG, *RHORMS;
	double *PAVG, *PRMS;
	double *TAVG, *TRMS;

	Options::StatsAvgType statsAvgType; 
	int stats_interval; 
	double prev_time;

	Stats(CurvilinearCSolver *cs, Options *opt){

	    this->cs  = cs;
	    this->opt = opt;
 
	    //get the flags for the statistics
	    velocityStats = opt->velocityStats;
	    thermoStats   = opt->thermoStats;

	    //If from restart file
	    if(opt->velocityStats){

		//allocate the arrays...

	        if(opt->velocityStatsFromRestart){
		     //read in the weight and the stat fields from the files

	        }else{
		
		     //reset the weights = 0.0
		     //initialize the fields
		}
	    }	


	    if(opt->thermoStats){

		//allocate the arrays

		//allocate the arrays...
		if(opt->thermoStatsFromRestart){
		    //read in the weight and the stat fields from the files

		}else{
		
		     //reset the weights = 0.0
		     //initialize the fields
		}

	    }


	    
	    if(opt->restartStats){
		weight = 0.0;
	    }
	}


};

#endif
