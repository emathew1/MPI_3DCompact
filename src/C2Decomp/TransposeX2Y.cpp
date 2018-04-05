#include "C2Decomp.hpp"

void C2Decomp::transposeX2Y(double *src, double *dst){

    int s1, s2, s3, d1, d2, d3;

    s1 = decompMain.xsz[0];
    s2 = decompMain.xsz[1];
    s3 = decompMain.xsz[2];

    d1 = decompMain.ysz[0];
    d2 = decompMain.ysz[1];
    d3 = decompMain.ysz[2];

    memSplitXY(src, s1, s2, s3, work1_r, dims[0], decompMain.x1dist);

    MPI_Alltoallv(work1_r, decompMain.x1cnts, decompMain.x1disp, realType, 
		  work2_r, decompMain.y1cnts, decompMain.y1disp, realType, 		   
		  DECOMP_2D_COMM_COL);

    memMergeXY(work2_r, d1, d2, d3, dst, dims[0], decompMain.y1dist);
    //memMergeXY_YMajor(work2_r, d1, d2, d3, dst, dims[0], decompMain.y1dist);

}


void C2Decomp::transposeX2Y_MajorIndex(double *src, double *dst){

    int s1, s2, s3, d1, d2, d3;

    s1 = decompMain.xsz[0];
    s2 = decompMain.xsz[1];
    s3 = decompMain.xsz[2];

    d1 = decompMain.ysz[0];
    d2 = decompMain.ysz[1];
    d3 = decompMain.ysz[2];

    //Always comes in major order...
    memSplitXY(src, s1, s2, s3, work1_r, dims[0], decompMain.x1dist);

    MPI_Alltoallv(work1_r, decompMain.x1cnts, decompMain.x1disp, realType, 
		  work2_r, decompMain.y1cnts, decompMain.y1disp, realType, 		   
		  DECOMP_2D_COMM_COL);

    memMergeXY_YMajor(work2_r, d1, d2, d3, dst, dims[0], decompMain.y1dist);

}

void C2Decomp::transposeX2Y_Start(MPI_Request &handle, double *src, double *dst, double *sbuf, double *rbuf){

    int s1, s2, s3;
   
    s1 = decompMain.xsz[0];
    s2 = decompMain.xsz[1];
    s3 = decompMain.xsz[2];

    memSplitXY(src, s1, s2, s3, sbuf, dims[0], decompMain.x1dist);
  
    MPI_Ialltoallv(sbuf, decompMain.x1cnts, decompMain.x1disp, realType, 
		   rbuf, decompMain.y1cnts, decompMain.y1disp, realType,
		   DECOMP_2D_COMM_COL, &handle);


} 

void C2Decomp::transposeX2Y_MajorIndex_Start(MPI_Request &handle, double *src, double *dst, double *sbuf, double *rbuf){

    int s1, s2, s3;
   
    s1 = decompMain.xsz[0];
    s2 = decompMain.xsz[1];
    s3 = decompMain.xsz[2];

    memSplitXY(src, s1, s2, s3, sbuf, dims[0], decompMain.x1dist);
  
    MPI_Ialltoallv(sbuf, decompMain.x1cnts, decompMain.x1disp, realType, 
		   rbuf, decompMain.y1cnts, decompMain.y1disp, realType,
		   DECOMP_2D_COMM_COL, &handle);


}

void C2Decomp::transposeX2Y_Wait(MPI_Request &handle, double *src, double *dst, double *sbuf, double *rbuf){

    int d1, d2, d3;
    MPI_Status status;

    d1 = decompMain.ysz[0];
    d2 = decompMain.ysz[1];
    d3 = decompMain.ysz[2];

    MPI_Wait(&handle, &status);

    memMergeXY(rbuf, d1, d2, d3, dst, dims[0], decompMain.y1dist);

}

void C2Decomp::transposeX2Y_MajorIndex_Wait(MPI_Request &handle, double *src, double *dst, double *sbuf, double *rbuf){

    int d1, d2, d3;
    MPI_Status status;

    d1 = decompMain.ysz[0];
    d2 = decompMain.ysz[1];
    d3 = decompMain.ysz[2];

    MPI_Wait(&handle, &status);

    memMergeXY_YMajor(rbuf, d1, d2, d3, dst, dims[0], decompMain.y1dist);

}
