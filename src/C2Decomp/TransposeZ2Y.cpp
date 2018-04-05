#include "C2Decomp.hpp"

void C2Decomp::transposeZ2Y(double *src, double *dst){

    int s1, s2, s3, d1, d2, d3;

    d1 = decompMain.ysz[0];
    d2 = decompMain.ysz[1];
    d3 = decompMain.ysz[2];


    MPI_Alltoallv(src,     decompMain.z2cnts, decompMain.z2disp, realType, 
		  work2_r, decompMain.y2cnts, decompMain.y2disp, realType, 
		  DECOMP_2D_COMM_ROW);

    
    memMergeZY(work2_r, d1, d2, d3, dst, dims[1], decompMain.y2dist);

}


void C2Decomp::transposeZ2Y_MajorIndex(double *src, double *dst){

    int s1, s2, s3, d1, d2, d3;

    s1 = decompMain.zsz[0];
    s2 = decompMain.zsz[1];
    s3 = decompMain.zsz[2];

    d1 = decompMain.ysz[0];
    d2 = decompMain.ysz[1];
    d3 = decompMain.ysz[2];

    //Just do the transpose here...
    for(int kp = 0; kp < s3; kp++){
        for(int jp = 0; jp < s2; jp++){
            for(int ip = 0; ip < s1; ip++){
                int ii  = kp*s2*s1 + jp*s1 + ip;
                int iip = jp*s3*s1 + ip*s3 + kp;
                work1_r[ii] = src[iip];
            }
        }
    }

    MPI_Alltoallv(work1_r, decompMain.z2cnts, decompMain.z2disp, realType, 
		  work2_r, decompMain.y2cnts, decompMain.y2disp, realType, 
		  DECOMP_2D_COMM_ROW);

    
    memMergeZY_YMajor(work2_r, d1, d2, d3, dst, dims[1], decompMain.y2dist);

}

void C2Decomp::transposeZ2Y_Start(MPI_Request &handle, double *src, double *dst, double *sbuf, double *rbuf){

    int s1, s2, s3;

    s1 = decompMain.zsz[0];
    s2 = decompMain.zsz[1];
    s3 = decompMain.zsz[2];

    memcpy(sbuf, src, s1*s2*s3*sizeof(double));

    MPI_Ialltoallv(sbuf, decompMain.z2cnts, decompMain.z2disp, realType,
                   rbuf, decompMain.y2cnts, decompMain.y2disp, realType,
                   DECOMP_2D_COMM_ROW, &handle);


}

void C2Decomp::transposeZ2Y_Wait(MPI_Request &handle, double *src, double *dst, double *sbuf, double *rbuf){

    int d1, d2, d3;
    MPI_Status status;

    d1 = decompMain.ysz[0];
    d2 = decompMain.ysz[1];
    d3 = decompMain.ysz[2];

    MPI_Wait(&handle, &status);

    memMergeZY(rbuf, d1, d2, d3, dst, dims[1], decompMain.y2dist);

}

void C2Decomp::transposeZ2Y_MajorIndex_Start(MPI_Request &handle, double *src, double *dst, double *sbuf, double *rbuf){

    int s1, s2, s3;

    s1 = decompMain.zsz[0];
    s2 = decompMain.zsz[1];
    s3 = decompMain.zsz[2];

    //Just do the transpose here...
    for(int kp = 0; kp < s3; kp++){
        for(int jp = 0; jp < s2; jp++){
            for(int ip = 0; ip < s1; ip++){
                int ii  = kp*s2*s1 + jp*s1 + ip;
                int iip = jp*s3*s1 + ip*s3 + kp;
                sbuf[ii] = src[iip];
            }
        }
    }

    MPI_Ialltoallv(sbuf, decompMain.z2cnts, decompMain.z2disp, realType,
                   rbuf, decompMain.y2cnts, decompMain.y2disp, realType,
                   DECOMP_2D_COMM_ROW, &handle);


}

void C2Decomp::transposeZ2Y_MajorIndex_Wait(MPI_Request &handle, double *src, double *dst, double *sbuf, double *rbuf){

    int d1, d2, d3;
    MPI_Status status;

    d1 = decompMain.ysz[0];
    d2 = decompMain.ysz[1];
    d3 = decompMain.ysz[2];

    MPI_Wait(&handle, &status);

    memMergeZY_YMajor(rbuf, d1, d2, d3, dst, dims[1], decompMain.y2dist);

}

