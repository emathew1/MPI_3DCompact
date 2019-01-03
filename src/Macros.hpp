#ifndef _MACROSH_
#define _MACROSH_

#define IF_RANK0 if(!mpiRank)

#define FOR_I3 for(int i = 0; i < 3; i++)
#define FOR_J3 for(int j = 0; j < 3; j++)
#define FOR_K3 for(int k = 0; k < 3; k++)

#define SIGNED_TET_VOLUME_6(A,B,C,D) (((B)[0]-(A)[0])*(((C)[1]-(A)[1])*((D)[2]-(A)[2])-((C)[2]-(A)[2])*((D)[1]-(A)[1])) + \
                                      ((B)[1]-(A)[1])*(((C)[2]-(A)[2])*((D)[0]-(A)[0])-((C)[0]-(A)[0])*((D)[2]-(A)[2])) + \
                                      ((B)[2]-(A)[2])*(((C)[0]-(A)[0])*((D)[1]-(A)[1])-((C)[1]-(A)[1])*((D)[0]-(A)[0])))

//Need to be sure that Nx,Ny,Nz are defined within scope to use these
//General indexing...
#define GETMAJIND_XPEN k*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETMAJIND_YPEN i*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_ZPEN j*pzSize[0]*pzSize[2] + i*pzSize[2] + k

#define GETIND_XPEN k*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETIND_YPEN k*pySize[0]*pySize[1] + j*pySize[0] + i
#define GETIND_ZPEN k*pzSize[0]*pzSize[1] + j*pzSize[0] + i

#define GETGLOBALXIND_XPEN pxStart[0] + i
#define GETGLOBALYIND_XPEN pxStart[1] + j
#define GETGLOBALZIND_XPEN pxStart[2] + k

#define GETGLOBALXIND_YPEN pyStart[0] + i
#define GETGLOBALYIND_YPEN pyStart[1] + j
#define GETGLOBALZIND_YPEN pyStart[2] + k

#define GETGLOBALXIND_ZPEN pzStart[0] + i
#define GETGLOBALYIND_ZPEN pzStart[1] + j
#define GETGLOBALZIND_ZPEN pzStart[2] + k


//For pulling data for Neumann BC's
#define GETMAJIND_XPEN_Xp1 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i+1
#define GETMAJIND_XPEN_Xp2 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i+2
#define GETMAJIND_XPEN_Xp3 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i+3
#define GETMAJIND_XPEN_Xp4 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i+4
#define GETMAJIND_XPEN_Xp5 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i+5
#define GETMAJIND_XPEN_Xp6 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i+6
#define GETMAJIND_XPEN_Xp7 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i+7
#define GETMAJIND_XPEN_Xp8 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i+8
#define GETMAJIND_XPEN_Xp9 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i+9
#define GETMAJIND_XPEN_Xp10 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i+10

#define GETMAJIND_XPEN_Xm1 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i-1
#define GETMAJIND_XPEN_Xm2 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i-2
#define GETMAJIND_XPEN_Xm3 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i-3
#define GETMAJIND_XPEN_Xm4 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i-4
#define GETMAJIND_XPEN_Xm5 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i-5
#define GETMAJIND_XPEN_Xm6 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i-6
#define GETMAJIND_XPEN_Xm7 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i-7
#define GETMAJIND_XPEN_Xm8 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i-8
#define GETMAJIND_XPEN_Xm9 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i-9
#define GETMAJIND_XPEN_Xm10 k*pxSize[0]*pxSize[1] + j*pxSize[0] + i-10

#define GETMAJIND_YPEN_Yp1 i*pySize[2]*pySize[1] + k*pySize[1] + j+1
#define GETMAJIND_YPEN_Yp2 i*pySize[2]*pySize[1] + k*pySize[1] + j+2
#define GETMAJIND_YPEN_Yp3 i*pySize[2]*pySize[1] + k*pySize[1] + j+3
#define GETMAJIND_YPEN_Yp4 i*pySize[2]*pySize[1] + k*pySize[1] + j+4
#define GETMAJIND_YPEN_Yp5 i*pySize[2]*pySize[1] + k*pySize[1] + j+5
#define GETMAJIND_YPEN_Yp6 i*pySize[2]*pySize[1] + k*pySize[1] + j+6
#define GETMAJIND_YPEN_Yp7 i*pySize[2]*pySize[1] + k*pySize[1] + j+7
#define GETMAJIND_YPEN_Yp8 i*pySize[2]*pySize[1] + k*pySize[1] + j+8
#define GETMAJIND_YPEN_Yp9 i*pySize[2]*pySize[1] + k*pySize[1] + j+9
#define GETMAJIND_YPEN_Yp10 i*pySize[2]*pySize[1] + k*pySize[1] + j+10

#define GETMAJIND_YPEN_Ym1 i*pySize[2]*pySize[1] + k*pySize[1] + j-1
#define GETMAJIND_YPEN_Ym2 i*pySize[2]*pySize[1] + k*pySize[1] + j-2
#define GETMAJIND_YPEN_Ym3 i*pySize[2]*pySize[1] + k*pySize[1] + j-3
#define GETMAJIND_YPEN_Ym4 i*pySize[2]*pySize[1] + k*pySize[1] + j-4
#define GETMAJIND_YPEN_Ym5 i*pySize[2]*pySize[1] + k*pySize[1] + j-5
#define GETMAJIND_YPEN_Ym6 i*pySize[2]*pySize[1] + k*pySize[1] + j-6
#define GETMAJIND_YPEN_Ym7 i*pySize[2]*pySize[1] + k*pySize[1] + j-7
#define GETMAJIND_YPEN_Ym8 i*pySize[2]*pySize[1] + k*pySize[1] + j-8
#define GETMAJIND_YPEN_Ym9 i*pySize[2]*pySize[1] + k*pySize[1] + j-9
#define GETMAJIND_YPEN_Ym10 i*pySize[2]*pySize[1] + k*pySize[1] + j-10

#define GETMAJIND_ZPEN_Zp1 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k+1
#define GETMAJIND_ZPEN_Zp2 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k+2
#define GETMAJIND_ZPEN_Zp3 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k+3
#define GETMAJIND_ZPEN_Zp4 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k+4
#define GETMAJIND_ZPEN_Zp5 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k+5
#define GETMAJIND_ZPEN_Zp6 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k+6
#define GETMAJIND_ZPEN_Zp7 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k+7
#define GETMAJIND_ZPEN_Zp8 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k+8
#define GETMAJIND_ZPEN_Zp9 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k+9
#define GETMAJIND_ZPEN_Zp10 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k+10

#define GETMAJIND_ZPEN_Zm1 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k-1
#define GETMAJIND_ZPEN_Zm2 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k-2
#define GETMAJIND_ZPEN_Zm3 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k-3
#define GETMAJIND_ZPEN_Zm4 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k-4
#define GETMAJIND_ZPEN_Zm5 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k-5
#define GETMAJIND_ZPEN_Zm6 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k-6
#define GETMAJIND_ZPEN_Zm7 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k-7
#define GETMAJIND_ZPEN_Zm8 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k-8
#define GETMAJIND_ZPEN_Zm9 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k-9
#define GETMAJIND_ZPEN_Zm10 j*pzSize[2]*pzSize[0] + i*pzSize[2] + k-10

#define GETYPIND_XPEN_Yp1 k*pxSize[0]*pxSize[1] + (j+1)*pxSize[0] + i
#define GETYPIND_XPEN_Yp2 k*pxSize[0]*pxSize[1] + (j+2)*pxSize[0] + i
#define GETYPIND_XPEN_Yp3 k*pxSize[0]*pxSize[1] + (j+3)*pxSize[0] + i
#define GETYPIND_XPEN_Yp4 k*pxSize[0]*pxSize[1] + (j+4)*pxSize[0] + i
#define GETYPIND_XPEN_Yp5 k*pxSize[0]*pxSize[1] + (j+5)*pxSize[0] + i
#define GETYPIND_XPEN_Yp6 k*pxSize[0]*pxSize[1] + (j+6)*pxSize[0] + i
#define GETYPIND_XPEN_Ym1 k*pxSize[0]*pxSize[1] + (j-1)*pxSize[0] + i
#define GETYPIND_XPEN_Ym2 k*pxSize[0]*pxSize[1] + (j-2)*pxSize[0] + i
#define GETYPIND_XPEN_Ym3 k*pxSize[0]*pxSize[1] + (j-3)*pxSize[0] + i
#define GETYPIND_XPEN_Ym4 k*pxSize[0]*pxSize[1] + (j-4)*pxSize[0] + i
#define GETYPIND_XPEN_Ym5 k*pxSize[0]*pxSize[1] + (j-5)*pxSize[0] + i
#define GETYPIND_XPEN_Ym6 k*pxSize[0]*pxSize[1] + (j-6)*pxSize[0] + i

#define GETZPIND_XPEN_Zp1 (k+1)*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETZPIND_XPEN_Zp2 (k+2)*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETZPIND_XPEN_Zp3 (k+3)*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETZPIND_XPEN_Zp4 (k+4)*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETZPIND_XPEN_Zp5 (k+5)*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETZPIND_XPEN_Zp6 (k+6)*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETZPIND_XPEN_Zm1 (k-1)*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETZPIND_XPEN_Zm2 (k-2)*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETZPIND_XPEN_Zm3 (k-3)*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETZPIND_XPEN_Zm4 (k-4)*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETZPIND_XPEN_Zm5 (k-5)*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETZPIND_XPEN_Zm6 (k-6)*pxSize[0]*pxSize[1] + j*pxSize[0] + i

#define GETMAJIND_YPEN_Xp1 (i+1)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xp2 (i+2)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xp3 (i+3)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xp4 (i+4)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xp5 (i+5)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xp6 (i+6)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xp7 (i+7)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xp8 (i+8)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xp9 (i+9)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xp10 (i+10)*pySize[2]*pySize[1] + k*pySize[1] + j

#define GETMAJIND_YPEN_Xm1 (i-1)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xm2 (i-2)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xm3 (i-3)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xm4 (i-4)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xm5 (i-5)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xm6 (i-6)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xm7 (i-7)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xm8 (i-8)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xm9 (i-9)*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_YPEN_Xm10 (i-10)*pySize[2]*pySize[1] + k*pySize[1] + j

#define GETMAJIND_YPEN_Zp1 i*pySize[2]*pySize[1] + (k+1)*pySize[1] + j
#define GETMAJIND_YPEN_Zp2 i*pySize[2]*pySize[1] + (k+2)*pySize[1] + j
#define GETMAJIND_YPEN_Zp3 i*pySize[2]*pySize[1] + (k+3)*pySize[1] + j
#define GETMAJIND_YPEN_Zp4 i*pySize[2]*pySize[1] + (k+4)*pySize[1] + j
#define GETMAJIND_YPEN_Zp5 i*pySize[2]*pySize[1] + (k+5)*pySize[1] + j
#define GETMAJIND_YPEN_Zp6 i*pySize[2]*pySize[1] + (k+6)*pySize[1] + j
#define GETMAJIND_YPEN_Zp7 i*pySize[2]*pySize[1] + (k+7)*pySize[1] + j
#define GETMAJIND_YPEN_Zp8 i*pySize[2]*pySize[1] + (k+8)*pySize[1] + j
#define GETMAJIND_YPEN_Zp9 i*pySize[2]*pySize[1] + (k+9)*pySize[1] + j
#define GETMAJIND_YPEN_Zp10 i*pySize[2]*pySize[1] + (k+10)*pySize[1] + j

#define GETMAJIND_YPEN_Zm1 i*pySize[2]*pySize[1] + (k-1)*pySize[1] + j
#define GETMAJIND_YPEN_Zm2 i*pySize[2]*pySize[1] + (k-2)*pySize[1] + j
#define GETMAJIND_YPEN_Zm3 i*pySize[2]*pySize[1] + (k-3)*pySize[1] + j
#define GETMAJIND_YPEN_Zm4 i*pySize[2]*pySize[1] + (k-4)*pySize[1] + j
#define GETMAJIND_YPEN_Zm5 i*pySize[2]*pySize[1] + (k-5)*pySize[1] + j
#define GETMAJIND_YPEN_Zm6 i*pySize[2]*pySize[1] + (k-6)*pySize[1] + j
#define GETMAJIND_YPEN_Zm7 i*pySize[2]*pySize[1] + (k-7)*pySize[1] + j
#define GETMAJIND_YPEN_Zm8 i*pySize[2]*pySize[1] + (k-8)*pySize[1] + j
#define GETMAJIND_YPEN_Zm9 i*pySize[2]*pySize[1] + (k-9)*pySize[1] + j
#define GETMAJIND_YPEN_Zm10 i*pySize[2]*pySize[1] + (k-10)*pySize[1] + j


#define FILL_GETMAJIND_YPEN_Xp {GETMAJIND_YPEN_Xp1, GETMAJIND_YPEN_Xp2, GETMAJIND_YPEN_Xp3, GETMAJIND_YPEN_Xp4, GETMAJIND_YPEN_Xp5, GETMAJIND_YPEN_Xp6, GETMAJIND_YPEN_Xp7, GETMAJIND_YPEN_Xp8, GETMAJIND_YPEN_Xp9, GETMAJIND_YPEN_Xp10}
#define FILL_GETMAJIND_YPEN_Xm {GETMAJIND_YPEN_Xm1, GETMAJIND_YPEN_Xm2, GETMAJIND_YPEN_Xm3, GETMAJIND_YPEN_Xm4, GETMAJIND_YPEN_Xm5, GETMAJIND_YPEN_Xm6, GETMAJIND_YPEN_Xm7, GETMAJIND_YPEN_Xm8, GETMAJIND_YPEN_Xm9, GETMAJIND_YPEN_Xm10}

#define FILL_GETMAJIND_YPEN_Yp {GETMAJIND_YPEN_Yp1, GETMAJIND_YPEN_Yp2, GETMAJIND_YPEN_Yp3, GETMAJIND_YPEN_Yp4, GETMAJIND_YPEN_Yp5, GETMAJIND_YPEN_Yp6, GETMAJIND_YPEN_Yp7, GETMAJIND_YPEN_Yp8, GETMAJIND_YPEN_Yp9, GETMAJIND_YPEN_Yp10}
#define FILL_GETMAJIND_YPEN_Ym {GETMAJIND_YPEN_Ym1, GETMAJIND_YPEN_Ym2, GETMAJIND_YPEN_Ym3, GETMAJIND_YPEN_Ym4, GETMAJIND_YPEN_Ym5, GETMAJIND_YPEN_Ym6, GETMAJIND_YPEN_Ym7, GETMAJIND_YPEN_Ym8, GETMAJIND_YPEN_Ym9, GETMAJIND_YPEN_Ym10}

#define FILL_GETMAJIND_YPEN_Zp {GETMAJIND_YPEN_Zp1, GETMAJIND_YPEN_Zp2, GETMAJIND_YPEN_Zp3, GETMAJIND_YPEN_Zp4, GETMAJIND_YPEN_Zp5, GETMAJIND_YPEN_Zp6, GETMAJIND_YPEN_Zp7, GETMAJIND_YPEN_Zp8, GETMAJIND_YPEN_Zp9, GETMAJIND_YPEN_Zp10}
#define FILL_GETMAJIND_YPEN_Zm {GETMAJIND_YPEN_Zm1, GETMAJIND_YPEN_Zm2, GETMAJIND_YPEN_Zm3, GETMAJIND_YPEN_Zm4, GETMAJIND_YPEN_Zm5, GETMAJIND_YPEN_Zm6, GETMAJIND_YPEN_Zm7, GETMAJIND_YPEN_Zm8, GETMAJIND_YPEN_Zm9, GETMAJIND_YPEN_Zm10}


//Loops
#define FOR_X_XPEN for(int i = 0; i < pxSize[0]; i++)
#define FOR_Y_XPEN for(int j = 0; j < pxSize[1]; j++)
#define FOR_Z_XPEN for(int k = 0; k < pxSize[2]; k++)
#define FOR_XYZ_XPEN for(int ip = 0; ip < pxSize[0]*pxSize[1]*pxSize[2]; ip++) 

#define FOR_X_YPEN for(int i = 0; i < pySize[0]; i++)
#define FOR_Y_YPEN for(int j = 0; j < pySize[1]; j++)
#define FOR_Z_YPEN for(int k = 0; k < pySize[2]; k++)
#define FOR_XYZ_YPEN for(int ip = 0; ip < pySize[0]*pySize[1]*pySize[2]; ip++) 

#define FOR_X_ZPEN for(int i = 0; i < pzSize[0]; i++)
#define FOR_Y_ZPEN for(int j = 0; j < pzSize[1]; j++)
#define FOR_Z_ZPEN for(int k = 0; k < pzSize[2]; k++)
#define FOR_XYZ_ZPEN for(int ip = 0; ip < pzSize[0]*pzSize[1]*pzSize[2]; ip++) 



//Need to close off each of these macros with their corresponding END_FOR macros'

//NOTE: These also need another check to see if the Start/End location are at global 0 or N-1
#define FOR_X0_XPEN {int i = 0;           if(GETGLOBALXIND_XPEN == 0){    for(int k = 0; k < pxSize[2]; k++){for(int j = 0; j < pxSize[1]; j++){int ip = GETIND_XPEN; 
#define FOR_X1_XPEN {int i = pxSize[0]-1; if(GETGLOBALXIND_XPEN == Nx-1){ for(int k = 0; k < pxSize[2]; k++){for(int j = 0; j < pxSize[1]; j++){int ip = GETIND_XPEN; 
#define FOR_Y0_XPEN {int j = 0;           if(GETGLOBALYIND_XPEN == 0){    for(int k = 0; k < pxSize[2]; k++){for(int i = 0; i < pxSize[0]; i++){int ip = GETIND_XPEN;
#define FOR_Y1_XPEN {int j = pxSize[1]-1; if(GETGLOBALYIND_XPEN == Ny-1){ for(int k = 0; k < pxSize[2]; k++){for(int i = 0; i < pxSize[0]; i++){int ip = GETIND_XPEN; 
#define FOR_Z0_XPEN {int k = 0;           if(GETGLOBALZIND_XPEN == 0){    for(int j = 0; j < pxSize[1]; j++){for(int i = 0; i < pxSize[0]; i++){int ip = GETIND_XPEN; 
#define FOR_Z1_XPEN {int k = pxSize[2]-1; if(GETGLOBALZIND_XPEN == Nz-1){ for(int j = 0; j < pxSize[1]; j++){for(int i = 0; i < pxSize[0]; i++){int ip = GETIND_XPEN; 

//NOTE: Major indexed are going to have the boundaries definitely included
#define FOR_X0_XPEN_MAJ FOR_X0_XPEN 
#define FOR_X1_XPEN_MAJ FOR_X1_XPEN
#define FOR_Y0_XPEN_MAJ FOR_Y0_XPEN
#define FOR_Y1_XPEN_MAJ FOR_Y1_XPEN
#define FOR_Z0_XPEN_MAJ FOR_Z0_XPEN
#define FOR_Z1_XPEN_MAJ FOR_Z1_XPEN

//NOTE: These also need another check to see if the Start/End location are at global 0 or N-1
#define FOR_X0_YPEN {int i = 0;           if(GETGLOBALXIND_YPEN == 0){    for(int k = 0; k < pySize[2]; k++){for(int j = 0; j < pySize[1]; j++){ int ip = GETIND_YPEN;  
#define FOR_X1_YPEN {int i = pySize[0]-1; if(GETGLOBALXIND_YPEN == Nx-1){ for(int k = 0; k < pySize[2]; k++){for(int j = 0; j < pySize[1]; j++){ int ip = GETIND_YPEN;  
#define FOR_Y0_YPEN {int j = 0;           if(GETGLOBALYIND_YPEN == 0){    for(int k = 0; k < pySize[2]; k++){for(int i = 0; i < pySize[0]; i++){ int ip = GETIND_YPEN;  
#define FOR_Y1_YPEN {int j = pySize[1]-1; if(GETGLOBALYIND_YPEN == Ny-1){ for(int k = 0; k < pySize[2]; k++){for(int i = 0; i < pySize[0]; i++){ int ip = GETIND_YPEN;  
#define FOR_Z0_YPEN {int k = 0;           if(GETGLOBALZIND_YPEN == 0){    for(int j = 0; j < pySize[1]; j++){for(int i = 0; i < pySize[0]; i++){ int ip = GETIND_YPEN;  
#define FOR_Z1_YPEN {int k = pySize[2]-1; if(GETGLOBALZIND_YPEN == Nz-1){ for(int j = 0; j < pySize[1]; j++){for(int i = 0; i < pySize[0]; i++){ int ip = GETIND_YPEN;  

//NOTE: Major indexed are going to have the boundaries definitely included
#define FOR_X0_YPEN_MAJ {int i = 0;           if(GETGLOBALXIND_YPEN == 0){    for(int k = 0; k < pySize[2]; k++){for(int j = 0; j < pySize[1]; j++){ int ip = GETMAJIND_YPEN; 
#define FOR_X1_YPEN_MAJ {int i = pySize[0]-1; if(GETGLOBALXIND_YPEN == Nx-1){ for(int k = 0; k < pySize[2]; k++){for(int j = 0; j < pySize[1]; j++){ int ip = GETMAJIND_YPEN;  
#define FOR_Y0_YPEN_MAJ {int j = 0;           if(GETGLOBALYIND_YPEN == 0){    for(int i = 0; i < pySize[0]; i++){for(int k = 0; k < pySize[2]; k++){ int ip = GETMAJIND_YPEN;  
#define FOR_Y1_YPEN_MAJ {int j = pySize[1]-1; if(GETGLOBALYIND_YPEN == Ny-1){ for(int i = 0; i < pySize[0]; i++){for(int k = 0; k < pySize[2]; k++){ int ip = GETMAJIND_YPEN;
#define FOR_Z0_YPEN_MAJ {int k = 0;           if(GETGLOBALZIND_YPEN == 0){    for(int i = 0; i < pySize[0]; i++){for(int j = 0; j < pySize[1]; j++){ int ip = GETMAJIND_YPEN; 
#define FOR_Z1_YPEN_MAJ {int k = pySize[2]-1; if(GETGLOBALZIND_YPEN == Nz-1){ for(int i = 0; i < pySize[0]; i++){for(int j = 0; j < pySize[1]; j++){ int ip = GETMAJIND_YPEN; 


//NOTE: These also need another check to see if the Start/End location are at global 0 or N-1
#define FOR_X0_ZPEN {int i = 0;           if(GETGLOBALXIND_ZPEN == 0){    for(int k = 0; k < pzSize[2]; k++){for(int j = 0; j < pzSize[1]; j++){int ip = GETIND_ZPEN;  
#define FOR_X1_ZPEN {int i = pzSize[0]-1; if(GETGLOBALXIND_ZPEN == Nx-1){ for(int k = 0; k < pzSize[2]; k++){for(int j = 0; j < pzSize[1]; j++){int ip = GETIND_ZPEN;
#define FOR_Y0_ZPEN {int j = 0;           if(GETGLOBALYIND_ZPEN == 0){    for(int k = 0; k < pzSize[2]; k++){for(int i = 0; i < pzSize[0]; i++){int ip = GETIND_ZPEN;
#define FOR_Y1_ZPEN {int j = pzSize[1]-1; if(GETGLOBALYIND_ZPEN == Ny-1){ for(int k = 0; k < pzSize[2]; k++){for(int i = 0; i < pzSize[0]; i++){int ip = GETIND_ZPEN;
#define FOR_Z0_ZPEN {int k = 0;           if(GETGLOBALZIND_ZPEN == 0){    for(int j = 0; j < pzSize[1]; j++){for(int i = 0; i < pzSize[0]; i++){int ip = GETIND_ZPEN; 
#define FOR_Z1_ZPEN {int k = pzSize[2]-1; if(GETGLOBALZIND_ZPEN == Nz-1){ for(int j = 0; j < pzSize[1]; j++){for(int i = 0; i < pzSize[0]; i++){int ip = GETIND_ZPEN;  

//NOTE: Major indexed are going to have the boundaries definitely included
#define FOR_X0_ZPEN_MAJ {int i = 0;           if(GETGLOBALXIND_ZPEN == 0){    for(int j = 0; j < pzSize[1]; j++){for(int k = 0; k < pzSize[2]; k++){int ip = GETMAJIND_ZPEN;  
#define FOR_X1_ZPEN_MAJ {int i = pzSize[0]-1; if(GETGLOBALXIND_ZPEN == Nx-1){ for(int j = 0; j < pzSize[1]; j++){for(int k = 0; k < pzSize[2]; k++){int ip = GETMAJIND_ZPEN;
#define FOR_Y0_ZPEN_MAJ {int j = 0;           if(GETGLOBALYIND_ZPEN == 0){    for(int i = 0; i < pzSize[0]; i++){for(int k = 0; k < pzSize[2]; k++){int ip = GETMAJIND_ZPEN;
#define FOR_Y1_ZPEN_MAJ {int j = pzSize[1]-1; if(GETGLOBALYIND_ZPEN == Ny-1){ for(int i = 0; i < pzSize[0]; i++){for(int k = 0; k < pzSize[2]; k++){int ip = GETMAJIND_ZPEN;
#define FOR_Z0_ZPEN_MAJ {int k = 0;           if(GETGLOBALZIND_ZPEN == 0){    for(int j = 0; j < pzSize[1]; j++){for(int i = 0; i < pzSize[0]; i++){int ip = GETMAJIND_ZPEN;
#define FOR_Z1_ZPEN_MAJ {int k = pzSize[2]-1; if(GETGLOBALZIND_ZPEN == Nz-1){ for(int j = 0; j < pzSize[1]; j++){for(int i = 0; i < pzSize[0]; i++){int ip = GETMAJIND_ZPEN;


#define END_FORX0 }}}}
#define END_FORX1 }}}}
#define END_FORY0 }}}}
#define END_FORY1 }}}}
#define END_FORZ0 }}}}
#define END_FORZ1 }}}}

#endif
