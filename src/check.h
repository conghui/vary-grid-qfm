#ifndef CHECK_H_NMIEIDTY
#define CHECK_H_NMIEIDTY

#include "fdutil.h"

float fsum1(const float *a, int n);
int isum1(const int *a, int n);
void write3df(const char *file, float ***dat, int n1, int n2, int n3);
void write2di(const char *file, int **dat, int n1, int n2);
void write2df(const char *file, float **dat, int n1, int n2);
void write1di(const char *file, int *dat, int n1);


void write3dfdm(const char *file, float ***dat, fdm3d fdm);

#endif /* end of include guard: CHECK_H_NMIEIDTY */
