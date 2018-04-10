#ifndef _MACROSH_
#define _MACROSH_

#define IF_RANK0 if(!mpiRank)

#define FOR_I3 for(int i = 0; i < 3; i++)
#define FOR_J3 for(int j = 0; j < 3; j++)
#define FOR_K3 for(int k = 0; k < 3; k++)

//Need to be sure that Nx,Ny,Nz are defined within scope to use these
//General indexing...
#define GETMAJIND_XPEN k*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETMAJIND_YPEN i*pySize[2]*pySize[1] + k*pySize[1] + j
#define GETMAJIND_ZPEN j*pzSize[0]*pzSize[2] + i*pzSize[2] + k

#define GETIND_XPEN k*pxSize[0]*pxSize[1] + j*pxSize[0] + i
#define GETIND_YPEN k*pySize[0]*pySize[1] + j*pySize[0] + i
#define GETIND_ZPEN k*pzSize[0]*pzSize[1] + j*pzSize[0] + i

//For pulling data for Neumann BC's
#define GET3DINDEX_XYZ_Xp1 k*Nx*Ny + j*Nx + i+1
#define GET3DINDEX_XYZ_Xp2 k*Nx*Ny + j*Nx + i+2
#define GET3DINDEX_XYZ_Xp3 k*Nx*Ny + j*Nx + i+3
#define GET3DINDEX_XYZ_Xp4 k*Nx*Ny + j*Nx + i+4
#define GET3DINDEX_XYZ_Xp5 k*Nx*Ny + j*Nx + i+5
#define GET3DINDEX_XYZ_Xp6 k*Nx*Ny + j*Nx + i+6
#define GET3DINDEX_XYZ_Xm1 k*Nx*Ny + j*Nx + i-1
#define GET3DINDEX_XYZ_Xm2 k*Nx*Ny + j*Nx + i-2
#define GET3DINDEX_XYZ_Xm3 k*Nx*Ny + j*Nx + i-3
#define GET3DINDEX_XYZ_Xm4 k*Nx*Ny + j*Nx + i-4
#define GET3DINDEX_XYZ_Xm5 k*Nx*Ny + j*Nx + i-5
#define GET3DINDEX_XYZ_Xm6 k*Nx*Ny + j*Nx + i-6

#define GET3DINDEX_XYZ_Yp1 k*Nx*Ny + (j+1)*Nx + i
#define GET3DINDEX_XYZ_Yp2 k*Nx*Ny + (j+2)*Nx + i
#define GET3DINDEX_XYZ_Yp3 k*Nx*Ny + (j+3)*Nx + i
#define GET3DINDEX_XYZ_Yp4 k*Nx*Ny + (j+4)*Nx + i
#define GET3DINDEX_XYZ_Yp5 k*Nx*Ny + (j+5)*Nx + i
#define GET3DINDEX_XYZ_Yp6 k*Nx*Ny + (j+6)*Nx + i
#define GET3DINDEX_XYZ_Ym1 k*Nx*Ny + (j-1)*Nx + i
#define GET3DINDEX_XYZ_Ym2 k*Nx*Ny + (j-2)*Nx + i
#define GET3DINDEX_XYZ_Ym3 k*Nx*Ny + (j-3)*Nx + i
#define GET3DINDEX_XYZ_Ym4 k*Nx*Ny + (j-4)*Nx + i
#define GET3DINDEX_XYZ_Ym5 k*Nx*Ny + (j-5)*Nx + i
#define GET3DINDEX_XYZ_Ym6 k*Nx*Ny + (j-6)*Nx + i

#define GET3DINDEX_XYZ_Zp1 (k+1)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zp2 (k+2)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zp3 (k+3)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zp4 (k+4)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zp5 (k+5)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zp6 (k+6)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zm1 (k-1)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zm2 (k-2)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zm3 (k-3)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zm4 (k-4)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zm5 (k-5)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zm6 (k-6)*Nx*Ny + j*Nx + i



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
#define FOR_X0_XPEN for(int k = 0; k < pxSize[2]; k++){for(int j = 0; j < pxSize[1]; j++){int i = 0; int ip = GETIND_XPEN; 
#define FOR_X1_XPEN for(int k = 0; k < pxSize[2]; k++){for(int j = 0; j < pxSize[1]; j++){int i = pxSize[0]-1; int ip = GETIND_XPEN; 
#define FOR_Y0_XPEN for(int k = 0; k < pxSize[2]; k++){for(int i = 0; i < pxSize[0]; i++){int j = 0; int ip = GETIND_XPEN; 
#define FOR_Y1_XPEN for(int k = 0; k < pxSize[2]; k++){for(int i = 0; i < pxSize[0]; i++){int j = pxSize[1]-1; int ip = GETIND_XPEN; 
#define FOR_Z0_XPEN for(int j = 0; j < pxSize[1]; j++){for(int i = 0; i < pxSize[0]; i++){int k = 0; int ip = GETIND_XPEN; 
#define FOR_Z1_XPEN for(int j = 0; j < pxSize[1]; j++){for(int i = 0; i < pxSize[0]; i++){int k = pxSize[2]-1; int ip = GETIND_XPEN; 

//NOTE: Major indexed are going to have the boundaries definitely included
#define FOR_X0_XPEN_MAJ FOR_X0_XPEN 
#define FOR_X1_XPEN_MAJ FOR_X0_XPEN
#define FOR_Y0_XPEN_MAJ FOR_Y0_XPEN
#define FOR_Y1_XPEN_MAJ FOR_Y1_XPEN
#define FOR_Z0_XPEN_MAJ FOR_Z0_XPEN
#define FOR_Z1_XPEN_MAJ FOR_Z1_XPEN

//NOTE: These also need another check to see if the Start/End location are at global 0 or N-1
#define FOR_X0_YPEN for(int k = 0; k < pySize[2]; k++){for(int j = 0; j < pySize[1]; j++){int i = 0; int ip = GETIND_YPEN; 
#define FOR_X1_YPEN for(int k = 0; k < pySize[2]; k++){for(int j = 0; j < pySize[1]; j++){int i = pySize[0]-1; int ip = GETIND_YPEN; 
#define FOR_Y0_YPEN for(int k = 0; k < pySize[2]; k++){for(int i = 0; i < pySize[0]; i++){int j = 0; int ip = GETIND_YPEN; 
#define FOR_Y1_YPEN for(int k = 0; k < pySize[2]; k++){for(int i = 0; i < pySize[0]; i++){int j = pySize[1]-1; int ip = GETIND_YPEN; 
#define FOR_Z0_YPEN for(int j = 0; j < pySize[1]; j++){for(int i = 0; i < pySize[0]; i++){int k = 0; int ip = GETIND_YPEN; 
#define FOR_Z1_YPEN for(int j = 0; j < pySize[1]; j++){for(int i = 0; i < pySize[0]; i++){int k = pySize[2]-1; int ip = GETIND_YPEN; 

//NOTE: Major indexed are going to have the boundaries definitely included
#define FOR_X0_YPEN_MAJ for(int k = 0; k < pySize[2]; k++){for(int j = 0; j < pySize[1]; j++){int i = 0; int ip = GETIND_YPEN_MAJ; 
#define FOR_X1_YPEN_MAJ for(int k = 0; k < pySize[2]; k++){for(int j = 0; j < pySize[1]; j++){int i = pySize[0]-1; int ip = GETIND_YPEN_MAJ; 
#define FOR_Y0_YPEN_MAJ for(int i = 0; i < pySize[0]; i++){for(int k = 0; k < pySize[2]; k++){int j = 0; int ip = GETIND_YPEN_MAJ; 
#define FOR_Y1_YPEN_MAJ for(int i = 0; i < pySize[0]; i++){for(int k = 0; k < pySize[2]; k++){int j = pySize[1]-1; int ip = GETIND_YPEN_MAJ; 
#define FOR_Z0_YPEN_MAJ for(int i = 0; i < pySize[0]; i++){for(int j = 0; j < pySize[1]; j++){int k = 0; int ip = GETIND_YPEN_MAJ; 
#define FOR_Z1_YPEN_MAJ for(int i = 0; i < pySize[0]; i++){for(int j = 0; j < pySize[1]; j++){int k = pySize[2]-1; int ip = GETIND_YPEN_MAJ; 


//NOTE: These also need another check to see if the Start/End location are at global 0 or N-1
#define FOR_X0_ZPEN for(int k = 0; k < pzSize[2]; k++){for(int j = 0; j < pzSize[1]; j++){int i = 0; int ip = GETIND_ZPEN; 
#define FOR_X1_ZPEN for(int k = 0; k < pzSize[2]; k++){for(int j = 0; j < pzSize[1]; j++){int i = pzSize[0]-1; int ip = GETIND_ZPEN; 
#define FOR_Y0_ZPEN for(int k = 0; k < pzSize[2]; k++){for(int i = 0; i < pzSize[0]; i++){int j = 0; int ip = GETIND_ZPEN; 
#define FOR_Y1_ZPEN for(int k = 0; k < pzSize[2]; k++){for(int i = 0; i < pzSize[0]; i++){int j = pzSize[1]-1; int ip = GETIND_ZPEN; 
#define FOR_Z0_ZPEN for(int j = 0; j < pzSize[1]; j++){for(int i = 0; i < pzSize[0]; i++){int k = 0; int ip = GETIND_ZPEN; 
#define FOR_Z1_ZPEN for(int j = 0; j < pzSize[1]; j++){for(int i = 0; i < pzSize[0]; i++){int k = pzSize[2]-1; int ip = GETIND_ZPEN; 

//NOTE: Major indexed are going to have the boundaries definitely included
#define FOR_X0_ZPEN_MAJ for(int j = 0; j < pzSize[1]; j++){for(int k = 0; k < pzSize[2]; k++){int i = 0; int ip = GETIND_ZPEN_MAJ; 
#define FOR_X1_ZPEN_MAJ for(int j = 0; j < pzSize[1]; j++){for(int k = 0; k < pzSize[2]; k++){int i = pzSize[0]-1; int ip = GETIND_ZPEN_MAJ; 
#define FOR_Y0_ZPEN_MAJ for(int i = 0; i < pzSize[0]; i++){for(int k = 0; k < pzSize[2]; k++){int j = 0; int ip = GETIND_ZPEN_MAJ; 
#define FOR_Y1_ZPEN_MAJ for(int i = 0; i < pzSize[0]; i++){for(int k = 0; k < pzSize[2]; k++){int j = pzSize[1]-1; int ip = GETIND_ZPEN_MAJ; 
#define FOR_Z0_ZPEN_MAJ for(int j = 0; j < pzSize[1]; j++){for(int i = 0; i < pzSize[0]; i++){int k = 0; int ip = GETIND_ZPEN_MAJ; 
#define FOR_Z1_ZPEN_MAJ for(int j = 0; j < pzSize[1]; j++){for(int i = 0; i < pzSize[0]; i++){int k = pzSize[2]-1; int ip = GETIND_ZPEN_MAJ; 



#define END_FORX0 }}
#define END_FORX1 }}
#define END_FORY0 }}
#define END_FORY1 }}
#define END_FORZ0 }}
#define END_FORZ1 }}

#endif
