#include "Stats.hpp"

void Stats::dumpStatsFields(){

};

void Stats::updateStatsFields(){

    double dt = cs->time-prev_time;
    
    if(velocityStats){

	double oldWeight     = velocityStatsWeight;
	double totNewWeight  = oldWeight + dt;

	FOR_XYZ_YPEN{
	    //Second order statistics...
	    URMS[ip]  = (oldWeight*(URMS[ip]*URMS[ip] + UAVG[ip]*UAVG[ip]) + dt*U[ip]*U[ip])/totNewWeight; 
	    VRMS[ip]  = (oldWeight*(VRMS[ip]*VRMS[ip] + VAVG[ip]*VAVG[ip]) + dt*V[ip]*V[ip])/totNewWeight; 
	    WRMS[ip]  = (oldWeight*(WRMS[ip]*WRMS[ip] + WAVG[ip]*WAVG[ip]) + dt*W[ip]*W[ip])/totNewWeight; 
	    UVREY[ip] = (oldWeight*(UVREY[ip]         + UAVG[ip]*VAVG[ip]) + dt*U[ip]*V[ip])/totNewWeight;
	    UWREY[ip] = (oldWeight*(UWREY[ip]         + UAVG[ip]*WAVG[ip]) + dt*U[ip]*W[ip])/totNewWeight;
	    VWREY[ip] = (oldWeight*(VWREY[ip]         + VAVG[ip]*WAVG[ip]) + dt*V[ip]*W[ip])/totNewWeight;

	    //First order stats
	    UAVG[ip]  = (oldWeight*UAVG[ip] + dt*U[ip])/totNewWeight;
	    VAVG[ip]  = (oldWeight*VAVG[ip] + dt*V[ip])/totNewWeight;
	    WAVG[ip]  = (oldWeight*WAVG[ip] + dt*W[ip])/totNewWeight;
	}

	//Do your averaging here...
	//doAveraging();

	//Finish up the second order stuff...
	FOR_XYZ_YPEN{	
	    URMS[ip] = sqrt(fabs(URMS[ip] - UAVG[ip]*UAVG[ip]));
	    VRMS[ip] = sqrt(fabs(VRMS[ip] - VAVG[ip]*VAVG[ip]));
	    WRMS[ip] = sqrt(fabs(WRMS[ip] - WAVG[ip]*WAVG[ip]));

	    UVREY[ip] = UVREY[ip] - UAVG[ip]*VAVG[ip];
	    UWREY[ip] = UWREY[ip] - UAVG[ip]*WAVG[ip];
	    VWREY[ip] = VWREY[ip] - VAVG[ip]*WAVG[ip];
	}	 	

    }


    if(thermoStats){

	double oldWeight     = thermoStatsWeight;
	double totNewWeight  = oldWeight + dt;

	FOR_XYZ_YPEN{
	    //Second order statistics...
	    RHORMS[ip]  = (oldWeight*(RHORMS[ip]*RHORMS[ip] + RHOAVG[ip]*RHOAVG[ip]) + dt*rho[ip]*rho[ip])/totNewWeight; 
	    PRMS[ip]    = (oldWeight*(PRMS[ip]*PRMS[ip] + PAVG[ip]*PAVG[ip]) + dt*p[ip]*p[ip])/totNewWeight; 
	    TRMS[ip]    = (oldWeight*(TRMS[ip]*TRMS[ip] + WAVG[ip]*TAVG[ip]) + dt*T[ip]*T[ip])/totNewWeight; 

	    //First order stats
	    RHOAVG[ip]  = (oldWeight*RHOAVG[ip] + dt*rho[ip])/totNewWeight;
	    PAVG[ip]    = (oldWeight*PAVG[ip]   + dt*p[ip])/totNewWeight;
	    TAVG[ip]    = (oldWeight*TAVG[ip]   + dt*T[ip])/totNewWeight;
	}

	//Do the averaging here...
	//doAveraging();

	//Finish up the second order stuff...
	FOR_XYZ_YPEN{	
	    RHORMS[ip] = sqrt(fabs(RHORMS[ip] - RHOAVG[ip]*RHOAVG[ip]));
	    PRMS[ip]   = sqrt(fabs(PRMS[ip]   - PAVG[ip]*PAVG[ip]));
	    TRMS[ip]   = sqrt(fabs(TRMS[ip]   - TAVG[ip]*TAVG[ip]));
	}	

    }


    prev_time = cs->time;
};

void Stats::doAveraging(double *phi){

    if(statsAvgType == Options::XI1_AVG){

	cs->c2d->transposeY2X_MajorIndex(phi, cs->tempX1);

	double *XI1_Avg = new double[pxSize[1]*pxSize[2]];

	for(int ip = 0; ip < pxSize[1]*pxSize[2]; ip++){
	    XI1_Avg[ip] = 0.0;
	}

	//Get the Xi index average
	FOR_Z_XPEN{
	    FOR_Y_XPEN{
		FOR_X_XPEN{

		    int ip = GETMAJIND_XPEN;
		    int jp = k*pxSize[1] + j;
 
	            XI1_Avg[jp] += cs->tempX1[ip]/(double)pxSize[0];
		}
	    }
	}

	//Load it back up into the container
	FOR_Z_XPEN{
	    FOR_Y_XPEN{
		FOR_X_XPEN{

		    int ip = GETMAJIND_XPEN;
		    int jp = k*pxSize[1] + j;
 		    
		    cs->tempX1[ip] = XI1_Avg[jp];
		}
	    }
	}

	delete[] XI1_Avg;

	cs->c2d->transposeX2Y_MajorIndex(cs->tempX1, phi);
				

    }else if(statsAvgType == Options::XI2_AVG){

	double *XI2_Avg = new double[pySize[0]*pySize[2]];

	for(int ip = 0; ip < pySize[0]*pySize[2]; ip++){
	    XI2_Avg[ip] = 0.0;
	}

	//Get the Xi index average
	FOR_Z_YPEN{
	    FOR_Y_YPEN{
		FOR_X_YPEN{

		    int ip = GETMAJIND_YPEN;
		    int jp = k*pySize[0] + i;
 
	            XI2_Avg[jp] += cs->tempY1[ip]/(double)pySize[1];
		}
	    }
	}

	//Load it back up into the container
	FOR_Z_YPEN{
	    FOR_Y_YPEN{
		FOR_X_YPEN{

		    int ip = GETMAJIND_YPEN;
		    int jp = k*pySize[0] + i;
 		    
		    cs->tempY1[ip] = XI2_Avg[jp];
		}
	    }
	}

	delete[] XI2_Avg;


    }else if(statsAvgType == Options::XI3_AVG){

	cs->c2d->transposeY2Z_MajorIndex(phi, cs->tempZ1);

	double *XI3_Avg = new double[pzSize[0]*pzSize[1]];

	for(int ip = 0; ip < pzSize[0]*pzSize[1]; ip++){
	    XI3_Avg[ip] = 0.0;
	}

	//Get the Xi index average
	FOR_Z_ZPEN{
	    FOR_Y_ZPEN{
		FOR_X_ZPEN{

		    int ip = GETMAJIND_ZPEN;
		    int jp = j*pzSize[0] + i;
 
	            XI3_Avg[jp] += cs->tempZ1[ip]/(double)pzSize[2];
		}
	    }
	}

	//Load it back up into the container
	FOR_Z_ZPEN{
	    FOR_Y_ZPEN{
		FOR_X_ZPEN{

		    int ip = GETMAJIND_ZPEN;
		    int jp = j*pzSize[0] + i;
 		    
		    cs->tempZ1[ip] = XI3_Avg[jp];
		}
	    }
	}

	delete[] XI3_Avg;

	cs->c2d->transposeZ2Y_MajorIndex(cs->tempZ1, phi);
	
    }else if(statsAvgType == Options::NONE){
	//just chill...
    }

}
