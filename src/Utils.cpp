#include "Utils.hpp"


void solveTri(double a[], double b[], double c[], double d[], double x[], double *work, int size)
{
        memcpy(work, b, size*sizeof(double));
        memcpy(x, d, size*sizeof(double));

        for(int ip = 1; ip < size; ip++){
            double m = a[ip]/work[ip-1];
            work[ip] = work[ip] - m*c[ip-1];
            x[ip] = x[ip] - m*x[ip-1];
        }

        x[size-1] /= work[size-1];
        for(int ip = size-2; ip >= 0; ip--){
            x[ip] = (x[ip] - c[ip]*x[ip+1])/work[ip];
        }

}

void cyclic(double *a, double *b, double *c, double alpha, double beta, double *r, int n, double *x)
{
        unsigned long i;
        double fact,gamma,bb[n],u[n],z[n],work[n];

        if (n <= 2) cout << "n too small in cyclic" << endl;

        gamma = -b[0];
        for (i=0;i<n;i++) bb[i]=b[i];
        bb[0]=b[0]-gamma;
        bb[n-1]=b[n-1]-alpha*beta/gamma;

        solveTri(a,bb,c,r,x,work,n);

        for (i=0;i<n;i++) u[i]=0.0;
        u[0]=gamma;
        u[n-1]=alpha;

        solveTri(a,bb,c,u,z,work,n);

        fact=(x[0]+beta*x[n-1]/gamma)/(1.0+z[0]+beta*z[n-1]/gamma);

        for (i=0;i<n;i++) x[i] -= fact*z[i];

}

void transposeMatrix(double *in, int Nx, int Ny, double *out){

    for(int i = 0; i < Ny; ++i)
        for(int j = 0; j < Nx; ++j)
            out[i*Nx + j] =  in[j*Ny + i];

}

void transposeMatrix_Fast1(const double *in, int n, int p, double *out, int block){
    for (int i = 0; i < n; i += block) {
        for(int j = 0; j < n; ++j) {
            for(int b = 0; b < block && i + b < n; ++b) {
                out[j*n + i + b] = in[(i + b)*n + j];
            }
        }
    }
}

void transposeMatrix_Fast2(const double *in, int n, int p, double *out, int blocksize){
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

void transposeXYZtoYZX(const double *in, int Nx, int Ny, int Nz, double *out){

    for(int ip = 0; ip < Nx; ip++){
	for(int kp = 0; kp < Nz; kp++){
	    for(int jp = 0; jp < Ny; jp++){
		out[ip*Nz*Ny + kp*Ny + jp] = in[kp*Ny*Nx + jp*Nx + ip];
	    }
	}
    }

} 

void transposeXYZtoYZX_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize){
    int i, j, k, row, col, sli;
    for ( i = 0; i < Nx; i += blocksize) {
	for( k = 0; k < Nz; k += blocksize){
            for ( j = 0; j < Ny; j += blocksize){
                for (row = i; row < i + blocksize && row < Nx; row++) {
		    for(sli = k; sli < k + blocksize && sli < Nz; sli++){
                    	for (col = j; col < j + blocksize && col < Ny; col++) {
                            out[row*Ny*Nz + sli*Ny + col] = in[sli*Ny*Nx + col*Nx + row];
			}
                    }
	        }
            }
        }
    }
}



void transposeYZXtoZXY(const double *in, int Nx, int Ny, int Nz, double *out){

    for(int jp = 0; jp < Ny; jp++){
    	for(int ip = 0; ip < Nx; ip++){
	    for(int kp = 0; kp < Nz; kp++){
		out[jp*Nx*Nz + ip*Nz + kp] = in[ip*Nz*Ny + kp*Ny + jp];
	    }
	}
    }

}

void transposeYZXtoZXY_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize){
    int i, j, k, row, col, sli;
    for ( j = 0; j < Ny; j += blocksize){
        for ( i = 0; i < Nx; i += blocksize) {
	    for( k = 0; k < Nz; k += blocksize){
                for (col = j; col < j + blocksize && col < Ny; col++) {
                    for (row = i; row < i + blocksize && row < Nx; row++) {
		        for(sli = k; sli < k + blocksize && sli < Nz; sli++){
                            out[col*Nx*Nz + row*Nz + sli] = in[row*Ny*Nz + sli*Ny + col];
			}
                    }
	        }
            }
        }
    }
}



void transposeXYZtoZXY(const double *in, int Nx, int Ny, int Nz, double *out){

    for(int jp = 0; jp < Ny; jp++){
    	for(int ip = 0; ip < Nx; ip++){
	    for(int kp = 0; kp < Nz; kp++){
		out[jp*Nx*Nz + ip*Nz + kp] = in[kp*Nx*Ny + jp*Nx + ip];
	    }
	}
    }

}


void transposeXYZtoZXY_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize){
    int i, j, k, row, col, sli;
    for ( j = 0; j < Ny; j += blocksize){
        for ( i = 0; i < Nx; i += blocksize) {
	    for( k = 0; k < Nz; k += blocksize){
                for (col = j; col < j + blocksize && col < Ny; col++) {
                    for (row = i; row < i + blocksize && row < Nx; row++) {
		        for(sli = k; sli < k + blocksize && sli < Nz; sli++){
                            out[col*Nx*Nz + row*Nz + sli] = in[sli*Ny*Nx + col*Nx + row];
			}
                    }
	        }
            }
        }
    }
}


void transposeZXYtoXYZ(const double *in, int Nx, int Ny, int Nz, double *out){

    for(int kp = 0; kp < Nz; kp++){
        for(int jp = 0; jp < Ny; jp++){
    	    for(int ip = 0; ip < Nx; ip++){
		out[kp*Nx*Ny + jp*Nx + ip] = in[jp*Nx*Nz + ip*Nz + kp];
	    }
	}
    }

}

void transposeZXYtoXYZ_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize){
    int i, j, k, row, col, sli;
    for( k = 0; k < Nz; k += blocksize){
        for ( j = 0; j < Ny; j += blocksize){
            for ( i = 0; i < Nx; i += blocksize) {
		for(sli = k; sli < k + blocksize && sli < Nz; sli++){
                    for (col = j; col < j + blocksize && col < Ny; col++) {
                        for (row = i; row < i + blocksize && row < Nx; row++) {
                            out[sli*Nx*Ny + col*Nx + row] = in[col*Nz*Nx + row*Nz + sli];
			}
                    }
	        }
            }
        }
    }
}

void transposeYZXtoXYZ(const double *in, int Nx, int Ny, int Nz, double *out){

    for(int kp = 0; kp < Nz; kp++){
    	for(int jp = 0; jp < Ny; jp++){
	    for(int ip = 0; ip < Nx; ip++){
		out[kp*Nx*Ny + jp*Nx + ip] = in[ip*Nz*Ny + kp*Ny + jp];
	    }
	}
    }
}

void transposeYZXtoXYZ_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize){
    int i, j, k, row, col, sli;
    for( k = 0; k < Nz; k += blocksize){
        for ( j = 0; j < Ny; j += blocksize){
            for ( i = 0; i < Nx; i += blocksize) {
		for(sli = k; sli < k + blocksize && sli < Nz; sli++){
                    for (col = j; col < j + blocksize && col < Ny; col++) {
                        for (row = i; row < i + blocksize && row < Nx; row++) {
                            out[sli*Nx*Ny + col*Nx + row] = in[row*Nz*Ny + sli*Ny + col];
			}
                    }
	        }
            }
        }
    }
}


double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void getRange(double *phi, string Var, int Nx, int Ny, int Nz, int mpiRank){

    double dataMin = 1000000.0;
    for(int ip = 0; ip < Nx*Ny*Nz; ip++)
	dataMin = min(phi[ip], dataMin);

    double dataMax = -1000000.0;
    for(int ip = 0; ip < Nx*Ny*Nz; ip++)
	dataMax = max(phi[ip], dataMax);

    double globalMin, globalMax;
    MPI_Reduce(&dataMin, &globalMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dataMax, &globalMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    IF_RANK0 cout << "  Range of " << Var << ": " << globalMin << ":" << globalMax << endl;

}

void getRangeValue(double *phi, int Nx, int Ny, int Nz, int mpiRank, double &globalDataMin, double &globalDataMax){

    double dataMin = 1000000.0;
    for(int ip = 0; ip < Nx*Ny*Nz; ip++)
	dataMin = min(phi[ip], dataMin);

    double dataMax = -1000000.0;
    for(int ip = 0; ip < Nx*Ny*Nz; ip++)
	dataMax = max(phi[ip], dataMax);

    MPI_Allreduce(&dataMin, &globalDataMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&dataMax, &globalDataMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);


}

bool isPointInHexa(double p[3], double vertex[8][3]){


    //Do bounding box check first
    double x_max = -10000000.0;
    double x_min =  10000000.0;
    double y_max = -10000000.0;
    double y_min =  10000000.0;
    double z_max = -10000000.0;
    double z_min =  10000000.0;

    for(int ip = 0; ip < 8; ip++){

	x_max = max(x_max, vertex[ip][0]);
	y_max = max(y_max, vertex[ip][1]);
	z_max = max(z_max, vertex[ip][2]);

	x_min = min(x_min, vertex[ip][0]);
	y_min = min(y_min, vertex[ip][1]);
	z_min = min(z_min, vertex[ip][2]);

    }

    bool isInBB = false;
    if(p[0] <= x_max && p[0] >= x_min){
	if(p[1] <= y_max && p[1] >= y_min){
	    if(p[2] <= z_max && p[2] >= z_min){
		isInBB = true;
	    }
	}
    }

    isInBB = true;
    //We are in the bounding box, continue with volume check
    if(isInBB){

	double volWithCentroid = getHexaVolume(vertex);
	double volWithPoint    = getHexaVolumeWithPoint(vertex, p);

	double eps = 1e-12;

	if(fabs(volWithPoint - volWithCentroid) < eps){
	    return true;
	}else{	
	    return false;
	}


    }else{
	return false;
    }

}

double getHexaVolume(double vertex[8][3]){

    //get hexa center...only guaranteed for convex polyhedra?
    double h_center[3] = {0,0,0};
    for(int ip = 0; ip < 8; ip++){
        FOR_I3 h_center[i] += vertex[ip][i];
    }
    FOR_I3 h_center[i] /= 8.0;

    //THESE ARE IN A SPECIFIC ORDER!!!
    double *A = vertex[0]; //{0.0, 0.0, 0.0} 
    double *B = vertex[2]; //{0.0, 1.0, 0.0}
    double *C = vertex[4]; //{1.0, 0.0, 0.0}
    double *D = vertex[6]; //{1.0, 1.0, 0.0}
    double *E = vertex[1]; //{0.0, 0.0, 1.0}
    double *F = vertex[3]; //{0.0, 1.0, 1.0}
    double *G = vertex[5]; //{1.0, 0.0, 1.0}
    double *H = vertex[7]; //{1.0, 1.0, 1.0}

    double volume_hexa = 0.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(A, B, C, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(B, D, C, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(G, F, E, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(F, G, H, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(E, B, A, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(B, E, F, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(C, D, G, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(H, G, D, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(A, C, G, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(C, G, E, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(F, D, B, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(D, F, H, h_center))/6.0;

    return volume_hexa;
};

double getHexaVolumeWithPoint(double vertex[8][3], double P[3]){

    //THESE ARE IN A SPECIFIC ORDER!!!
    double *A = vertex[0]; //{0.0, 0.0, 0.0} 
    double *B = vertex[2]; //{0.0, 1.0, 0.0}
    double *C = vertex[4]; //{1.0, 0.0, 0.0}
    double *D = vertex[6]; //{1.0, 1.0, 0.0}
    double *E = vertex[1]; //{0.0, 0.0, 1.0}
    double *F = vertex[3]; //{0.0, 1.0, 1.0}
    double *G = vertex[5]; //{1.0, 0.0, 1.0}
    double *H = vertex[7]; //{1.0, 1.0, 1.0}

    double volume_hexa = 0.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(A, B, C, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(B, D, C, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(G, F, E, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(F, G, H, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(E, B, A, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(B, E, F, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(C, D, G, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(H, G, D, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(A, C, G, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(C, G, E, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(F, D, B, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(D, F, H, P))/6.0;



    return volume_hexa;
};
