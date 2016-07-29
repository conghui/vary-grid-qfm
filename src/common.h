#ifndef COMMON_H_XWGHYQZP
#define COMMON_H_XWGHYQZP
#include <rsf.h>
#include "fdutil.h"

#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
#endif

 float* compute_fdcoef(int nop, float dz, float dx, float dy,
    bool is_optimized, const int d_order);
 int factorial(int n);
 float* normal_fdcoef(int nop, const int d_order);
 float* optimal_fdcoef(int nop, const int d_order);

 float* damp_make(int ntransit);
 void apply_abc(float*** uu2, float*** uu1, int nz, int nx, int ny, int nbd,
    abcone3d abc, int nop, float* damp);


#endif /* end of include guard: COMMON_H_XWGHYQZP */

