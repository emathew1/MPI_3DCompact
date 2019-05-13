#include "AbstractSGS.hpp"

//The default is the assume we're doing mu_sgs averaging with the default filters...
void AbstractSGS::doAveraging(){
    doAveraging(mu_sgs, mu_sgs, filtX, filtY, filtZ);
}

//Next is just to assume we're using the default filters...
void AbstractSGS::doAveraging(double *phi, double *phiF){
    doAveraging(phi, phiF, filtX, filtY, filtZ);
}

void AbstractSGS::doAveraging(double *phi, double *phiF, AbstractFilter *fx, AbstractFilter *fy, AbstractFilter *fz){

    if(musgsAvgType == Options::XI1_AVG){

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

	cs->c2d->transposeX2Y_MajorIndex(cs->tempX1, phiF);
				

    }else if(musgsAvgType == Options::XI2_AVG){

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
 
	            XI2_Avg[jp] += phi[ip]/(double)pySize[1];
		}
	    }
	}

	//Load it back up into the container
	FOR_Z_YPEN{
	    FOR_Y_YPEN{
		FOR_X_YPEN{

		    int ip = GETMAJIND_YPEN;
		    int jp = k*pySize[0] + i;
 		    
		    phiF[ip] = XI2_Avg[jp];
		}
	    }
	}

	delete[] XI2_Avg;


    }else if(musgsAvgType == Options::XI3_AVG){

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

	cs->c2d->transposeZ2Y_MajorIndex(cs->tempZ1, phiF);

    }else if(musgsAvgType == Options::GLOBAL){

        cout << "NEED TO IMPLEMENT GLOBAL MUSGS AVERAGING" << endl;

    }else if(musgsAvgType == Options::LOCAL){

	//Should just impletement a legit filter
	fy->filterField(phi, cs->tempY1);
	cs->c2d->transposeY2Z_MajorIndex(cs->tempY1, cs->tempZ1);
	fz->filterField(cs->tempZ1, cs->tempZ2);
	cs->c2d->transposeZ2Y_MajorIndex(cs->tempZ2, cs->tempY1);
	cs->c2d->transposeY2X_MajorIndex(cs->tempY1, cs->tempX1);
	fx->filterField(cs->tempX1, cs->tempX2);
	cs->c2d->transposeX2Y_MajorIndex(cs->tempX2, phiF);

    }else if(musgsAvgType == Options::NONE){
	//just chill...
    }

}
