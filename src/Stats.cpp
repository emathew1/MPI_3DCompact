#include "Stats.hpp"

void Stats::dumpStatsFields(){

    if(velocityStats){

        //Convert over to X-major indexing
	FOR_Z_YPEN{
	    FOR_Y_YPEN{
		FOR_X_YPEN{
		    int ip = GETMAJIND_YPEN;
		    int jp = GETIND_YPEN;

		    cs->tempY1[jp] = UAVG[ip];
		    cs->tempY2[jp] = VAVG[ip];
		    cs->tempY3[jp] = WAVG[ip];
		    cs->tempY4[jp] = URMS[ip];
		    cs->tempY5[jp] = VRMS[ip];
		    cs->tempY6[jp] = WRMS[ip];
		    cs->tempY7[jp] = UVREY[ip];
		    cs->tempY8[jp] = UWREY[ip];
		    cs->tempY9[jp] = VWREY[ip];
		}
	    }
	}

	IF_RANK0{
            cout << endl;
            cout << " > ========================" << endl;
            cout << " >  DUMPING VELOCITY STATS " << endl;
            cout << " > ========================" << endl;
	}

	ofstream outfile;
	outfile.precision(17);
	string outputFilename;
	outputFilename = "VelocityStats.";
	ostringstream timeStepString;
	timeStepString << cs->timeStep;

	outputFilename.append(timeStepString.str());

	MPI_File fh;
	MPI_Offset disp, filesize;

	MPI_File_open(MPI_COMM_WORLD, outputFilename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	filesize = 0;
	MPI_File_set_size(fh, filesize);

	disp = 0;

	cs->c2d->writeScalar(fh, disp, 1, &velocityStatsWeight);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY1);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY2);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY3);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY4);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY5);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY6);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY7);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY8);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY9);

	MPI_File_close(&fh);

    }

    if(thermoStats){

        //Convert over to X-major indexing
	FOR_Z_YPEN{
	    FOR_Y_YPEN{
		FOR_X_YPEN{
		    int ip = GETMAJIND_YPEN;
		    int jp = GETIND_YPEN;

		    cs->tempY1[jp] = RHOAVG[ip];
		    cs->tempY2[jp] = RHORMS[ip];
		    cs->tempY3[jp] = PAVG[ip];
		    cs->tempY4[jp] = PRMS[ip];
		    cs->tempY5[jp] = TAVG[ip];
		    cs->tempY6[jp] = TRMS[ip];
		}
	    }
	}

	IF_RANK0{
            cout << endl;
            cout << " > ======================" << endl;
            cout << " >  DUMPING THERMO STATS " << endl;
            cout << " > ======================" << endl;
	}

	ofstream outfile;
	outfile.precision(17);
	string outputFilename;
	outputFilename = "ThermoStats.";
	ostringstream timeStepString;
	timeStepString << cs->timeStep;

	outputFilename.append(timeStepString.str());

	MPI_File fh;
	MPI_Offset disp, filesize;

	MPI_File_open(MPI_COMM_WORLD, outputFilename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	filesize = 0;
	MPI_File_set_size(fh, filesize);

	disp = 0;

	cs->c2d->writeScalar(fh, disp, 1, &thermoStatsWeight);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY1);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY2);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY3);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY4);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY5);
	cs->c2d->writeVar(fh, disp, 1, cs->tempY6);

	MPI_File_close(&fh);

    }

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
	doAveraging(URMS);
	doAveraging(VRMS);
	doAveraging(WRMS);
	doAveraging(UVREY);
	doAveraging(UWREY);
	doAveraging(VWREY);
	doAveraging(UAVG);
	doAveraging(VAVG);
	doAveraging(WAVG);

	//Finish up the second order stuff...
	FOR_XYZ_YPEN{	
	    URMS[ip] = sqrt(fabs(URMS[ip] - UAVG[ip]*UAVG[ip]));
	    VRMS[ip] = sqrt(fabs(VRMS[ip] - VAVG[ip]*VAVG[ip]));
	    WRMS[ip] = sqrt(fabs(WRMS[ip] - WAVG[ip]*WAVG[ip]));

	    UVREY[ip] = UVREY[ip] - UAVG[ip]*VAVG[ip];
	    UWREY[ip] = UWREY[ip] - UAVG[ip]*WAVG[ip];
	    VWREY[ip] = VWREY[ip] - VAVG[ip]*WAVG[ip];
	}	 	

	IF_RANK0{
	    cout << " > Updated velocity stats, current stats weight = " << totNewWeight << endl;
	}

        velocityStatsWeight += dt;
    }


    if(thermoStats){

	double oldWeight     = thermoStatsWeight;
	double totNewWeight  = oldWeight + dt;

	FOR_XYZ_YPEN{
	    //Second order statistics...
	    RHORMS[ip]  = (oldWeight*(RHORMS[ip]*RHORMS[ip] + RHOAVG[ip]*RHOAVG[ip]) + dt*rho[ip]*rho[ip])/totNewWeight; 
	    PRMS[ip]    = (oldWeight*(PRMS[ip]*PRMS[ip] + PAVG[ip]*PAVG[ip]) + dt*p[ip]*p[ip])/totNewWeight; 
	    TRMS[ip]    = (oldWeight*(TRMS[ip]*TRMS[ip] + TAVG[ip]*TAVG[ip]) + dt*T[ip]*T[ip])/totNewWeight; 

	    //First order stats
	    RHOAVG[ip]  = (oldWeight*RHOAVG[ip] + dt*rho[ip])/totNewWeight;
	    PAVG[ip]    = (oldWeight*PAVG[ip]   + dt*p[ip])/totNewWeight;
	    TAVG[ip]    = (oldWeight*TAVG[ip]   + dt*T[ip])/totNewWeight;
	}

	//Do the averaging here...
	doAveraging(RHORMS);
	doAveraging(PRMS);
	doAveraging(TRMS);
	doAveraging(RHOAVG);
	doAveraging(PAVG);
	doAveraging(TAVG);

	//Finish up the second order stuff...
	FOR_XYZ_YPEN{	
	    RHORMS[ip] = sqrt(fabs(RHORMS[ip] - RHOAVG[ip]*RHOAVG[ip]));
	    PRMS[ip]   = sqrt(fabs(PRMS[ip]   - PAVG[ip]*PAVG[ip]));
	    TRMS[ip]   = sqrt(fabs(TRMS[ip]   - TAVG[ip]*TAVG[ip]));
	}	

	IF_RANK0{
	    cout << " > Updated thermo stats, current stats weight = " << totNewWeight << endl;
	}

        thermoStatsWeight += dt;
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
