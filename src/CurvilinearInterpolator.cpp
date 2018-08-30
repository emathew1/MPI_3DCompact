#include "CurvilinearInterpolator.hpp"

void CurvilinearInterpolator::getDataHalo(double *dataIn, double *dataHalo){

     //need to update halo data on this dataIn...
     double *data_temp1;
     cs->c2d->allocX(data_temp1);
     cs->c2d->transposeY2X_MajorIndex(dataIn, data_temp1);

     double *data_temp2;
     cs->c2d->allocY(data_temp2);
     cs->c2d->transposeX2Y(data_temp1, data_temp2);

     delete[] data_temp1;

     cs->c2d->updateHalo(data_temp2, data_halo, 1, 1);
     delete[] data_temp2;

};
