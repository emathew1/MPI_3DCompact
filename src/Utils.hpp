#ifndef _UTILSH_
#define _UTILSH_

#include <math.h>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "Macros.hpp"
#include "Domain.hpp"

using namespace std;

//Tri/Pentadiagonal Matrix solver stuff
void solveTri(double a[], double b[], double c[], double d[], double x[], double *work, int n);
void cyclic(double *a, double *b, double *c, double alpha, double beta,
        double *r, int n, double *x);
void solvePenta(double e[], double a[], double d[], double c[], double f[], double b[], double x[], double B[], double C[], double D[], int size); 
void cyclicPenta(double a[], double b[], double c[], double d[], double e[], double f[], double cp1[6], double x[], int N);

//Stuff to invert 4x4 matrix needed in cyclicPenta
double MINOR(double m[16], int r0, int r1, int r2, int c0, int c1, int c2);
void adjoint(double m[16], double adjOut[16]);
double det(double m[16]);
void invertRowMajor(double m[16], double invOut[16]);

//Transpose Stuff
void transposeMatrix(double *in, int Nx, int Ny, double *out);
void transposeMatrix_Fast1(const double *in, int n, int p, double *out, int block);
void transposeMatrix_Fast2(const double *in, int n, int p, double *out, int blocksize);

void transposeXYZtoYZX(const double *in, int Nx, int Ny, int Nz, double *out);
void transposeXYZtoYZX_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize);
void transposeXYZtoZXY(const double *in, int Nx, int Ny, int Nz, double *out);
void transposeXYZtoZXY_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize);
void transposeYZXtoZXY(const double *in, int Nx, int Ny, int Nz, double *out);
void transposeYZXtoZXY_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize);
void transposeZXYtoXYZ(const double *in, int Nx, int Ny, int Nz, double *out);
void transposeZXYtoXYZ_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize);
void transposeYZXtoXYZ(const double *in, int Nx, int Ny, int Nz, double *out);
void transposeYZXtoXYZ_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize);

//Random number generation
double fRand(double fMin, double fMax);


void getRange(double *phi, string Var, int Nx, int Ny, int Nz, int mpiRank);
void getRangeValue(double *phi, int Nx, int Ny, int Nz, int mpiRank, double &globalDataMin, double &globalDataMax);

bool isPointInHexa(double p[3], double vertex[8][3]);
double getHexaVolume(double vertex[8][3]);
double getHexaVolumeWithPoint(double vertex[8][3], double P[3]);

void getDataFromIndex(double *f_in, int *index, int N, double *f_out);

#endif
