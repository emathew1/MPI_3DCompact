#include "C2Decomp.hpp"

void C2Decomp::fastLocalTranspose(const double *in, int n, int p, double *out, int blocksize){
    int i, j, row, col;
    for ( i = 0; i < n; i += blocksize) {
        for ( j = 0; j < p; j += blocksize) {
            for (row = i; row < i + blocksize && row < n; row++) {
                for (col = j; col < j + blocksize && col < p; col++) {
                    out[row*p + col] = in[col*n + row];
                }
            }
        }
    }
}

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


void C2Decomp::transposeChunkX2Y_MajorIndex_Start(MPI_Request &handle, double **src, double **dst, double *sbuf, double *rbuf, int numArrayInChunk, double *src2){

    int s1, s2, s3;
   
    s1 = decompMain.xsz[0];
    s2 = decompMain.xsz[1];
    s3 = decompMain.xsz[2];

    MPI_Datatype MPI_DOUBLECHUNK;
    MPI_Type_contiguous(numArrayInChunk, MPI_DOUBLE, &MPI_DOUBLECHUNK);
    MPI_Type_commit(&MPI_DOUBLECHUNK);

    for(int ip = 0; ip < numArrayInChunk; ip++){
	memSplitXY(src[ip], s1, s2, s3, &src2[ip*s1*s2*s3], dims[0], decompMain.x1dist);
    }

    fastLocalTranspose(src2, s1*s2*s3, numArrayInChunk, sbuf, 8);
/*
    //Really just a transpose...
    for(int ip = 0; ip < s1*s2*s3; ip++){
	for(int jp = 0; jp < numArrayInChunk; jp++){
	    sbuf[ip*numArrayInChunk + jp] = src2[jp*s1*s2*s3 + ip];
	}
    } 
  */
    MPI_Ialltoallv(sbuf, decompMain.x1cnts, decompMain.x1disp, MPI_DOUBLECHUNK, 
		   rbuf, decompMain.y1cnts, decompMain.y1disp, MPI_DOUBLECHUNK,
		   DECOMP_2D_COMM_COL, &handle);


}

void C2Decomp::transposeChunkX2Y_MajorIndex_Wait(MPI_Request &handle, double **src, double **dst, double *sbuf, double *rbuf, int numArrayInChunk, double *dst2){

    int d1, d2, d3;
    MPI_Status status;

    d1 = decompMain.ysz[0];
    d2 = decompMain.ysz[1];
    d3 = decompMain.ysz[2];

    MPI_Wait(&handle, &status);


    fastLocalTranspose(rbuf, numArrayInChunk, d1*d2*d3, dst2, 8);
/*
    //Again just a transpose...
    for(int ip = 0; ip < d1*d2*d3; ip++){
	for(int jp = 0; jp < numArrayInChunk; jp++){
	    dst2[jp*d1*d2*d3 + ip] = rbuf[ip*numArrayInChunk + jp];
	}
    }
*/
    for(int ip = 0; ip < numArrayInChunk; ip++){
        memMergeXY_YMajor(&dst2[ip*d1*d2*d3], d1, d2, d3, dst[ip], dims[0], decompMain.y1dist);
    }

}


void C2Decomp::transposeChunkX2Y_MajorIndex(double **src, double **dst, int numArrayInChunk, double *src2, double *dst2){

    int s1, s2, s3;
    int d1, d2, d3;  
 
    s1 = decompMain.xsz[0];
    s2 = decompMain.xsz[1];
    s3 = decompMain.xsz[2];
    d1 = decompMain.ysz[0];
    d2 = decompMain.ysz[1];
    d3 = decompMain.ysz[2];

    MPI_Datatype MPI_DOUBLECHUNK;
    MPI_Type_contiguous(numArrayInChunk, MPI_DOUBLE, &MPI_DOUBLECHUNK);
    MPI_Type_commit(&MPI_DOUBLECHUNK);

    for(int ip = 0; ip < numArrayInChunk; ip++){
	memSplitXY(src[ip], s1, s2, s3, &src2[ip*s1*s2*s3], dims[0], decompMain.x1dist);
    }

    double *work1 = new double[numArrayInChunk*s1*s2*s3];
    double *work2 = new double[numArrayInChunk*d1*d2*d3];

    fastLocalTranspose(src2, s1*s2*s3, numArrayInChunk, work1, 8);
    
    MPI_Alltoallv(work1, decompMain.x1cnts, decompMain.x1disp, MPI_DOUBLECHUNK, 
		  work2, decompMain.y1cnts, decompMain.y1disp, MPI_DOUBLECHUNK,
		  DECOMP_2D_COMM_COL);

    fastLocalTranspose(work2, numArrayInChunk, d1*d2*d3, dst2, 8);
 
    for(int ip = 0; ip < numArrayInChunk; ip++){
        memMergeXY_YMajor(&dst2[ip*d1*d2*d3], d1, d2, d3, dst[ip], dims[0], decompMain.y1dist);
    }

}



