#include "CurvilinearInterpolator.hpp"

void CurvilinearInterpolator::getDataHalo(double *dataIn, double *&dataHalo){

     //need to update halo data on this dataIn...
     double *data_temp1;
     cs->c2d->allocX(data_temp1);
     cs->c2d->transposeY2X_MajorIndex(dataIn, data_temp1);

     double *data_temp2;
     cs->c2d->allocY(data_temp2);
     cs->c2d->transposeX2Y(data_temp1, data_temp2);

     delete[] data_temp1;

     cs->c2d->updateHalo(data_temp2, dataHalo, 1, 1);
     delete[] data_temp2;

};

void CurvilinearInterpolator::getOrderedBlockData(int ip, int jp, int kp, double *dataHalo, double box_data[8]){


     int iih_0_0_0;
     int iih_0_0_1;
     int iih_0_1_0;
     int iih_0_1_1;
     int iih_1_0_0;
     int iih_1_0_1;
     int iih_1_1_0;
     int iih_1_1_1;

     //This is the halo array index for the same point
     iih_0_0_0 = (kp+1)*pySize[1]*(pySize[0]+2) + jp*(pySize[0]+2) + ip + 1;

     //Halo array index for i, j, k+1
     iih_0_0_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (jp)*(pySize[0]+2) + ip + 1;

     //Halo array index for i, j+1, k
     if(transPeriodicY && jp == (Ny-1)){
         iih_0_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (0)*(pySize[0]+2) + ip + 1;
     }else{
         iih_0_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 1;
     }

     //Halo array index for i, j+1, k+1
     if(transPeriodicY && jp == (Ny-1)){
         iih_0_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (0)*(pySize[0]+2) + ip + 1;
     }else{
         iih_0_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 1;
     }

     //Halo array index for i+1, j, k
     iih_1_0_0 = (kp+1)*pySize[1]*(pySize[0]+2) + jp*(pySize[0]+2) + ip + 2;

     //Halo array index for i+1, j, k+1
     iih_1_0_1 = (kp+2)*pySize[1]*(pySize[0]+2) + jp*(pySize[0]+2) + ip + 2;

     //Halo array index for i+1, j+1, k
     if(transPeriodicY && jp == (Ny-1)){
         iih_1_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (0)*(pySize[0]+2) + ip + 2;
     }else{
         iih_1_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 2;
     }

     //Halo array index for i+1, j+1, k+1
     if(transPeriodicY && jp == (Ny-1)){
         iih_1_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (0)*(pySize[0]+2) + ip + 2;
     }else{
         iih_1_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 2;
     }


     box_data[0] = dataHalo[iih_0_0_0];
     box_data[1] = dataHalo[iih_0_0_1];
     box_data[2] = dataHalo[iih_0_1_0]; 
     box_data[3] = dataHalo[iih_0_1_1];
     box_data[4] = dataHalo[iih_1_0_0];
     box_data[5] = dataHalo[iih_1_0_1];
     box_data[6] = dataHalo[iih_1_1_0];
     box_data[7] = dataHalo[iih_1_1_1];

};

void CurvilinearInterpolator::interpolateData(double *dataIn, double *interpedDataOut){

    //Get the halo of the data
    double *dataHalo = NULL;
    getDataHalo(dataIn, dataHalo);

    for(int ii = 0; ii < localPointFoundCount; ii++){

        //Get the ip, jp, kp based off index
        int icv = icvList[ii];
        int jp = icv%pySize[1];
        int kp = (icv/pySize[1])%pySize[2];
        int ip =  icv/(pySize[2]*pySize[1]);

	//Pull the required data from the halo array
	double box_data[8];
	getOrderedBlockData(ip, jp, kp, dataHalo, box_data);

	//Interpolate using the weights for this point...
	interpedDataOut[ii] = 0;
	for(int jj = 0; jj < 8; jj++){
	    interpedDataOut[ii] += Ni[ii][jj]*box_data[jj]; 
	}
    }
  
    delete[] dataHalo;
};
