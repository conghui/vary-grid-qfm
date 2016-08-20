#include "common.h"



#ifndef SF_HAS_SSE
#undef __SSE__
#endif
#ifndef SF_HAS_AVX
#undef __AVX__
#endif

#ifdef sun
#define restrict
#endif

#if defined __SSE__ || defined __AVX__
#include <immintrin.h>
#endif


  float*
    damp_make(int ntransit)
    {
      int ib;
      float* damp = NULL;
      float sb, fb;
      if (ntransit>0) damp = sf_floatalloc(ntransit);
      sb = 4.0*ntransit;
      for(ib=0; ib<ntransit; ib++) {
        fb = ib/(sqrt(2.0)*sb);
        damp[ib] = exp(-fb*fb);
      }
      return damp;
    }

  float*
    compute_fdcoef(int nop, float dz, float dx, float dy, bool is_optimized, const int d_order)
    /*
       Optimized fd coeffifients from:
       1. Yang Liu. "Globally optimal finite-difference schemes based on least squares" Geophysics 78.4 (2013)
       2. Zhang, Jin-Hai, and Zhen-Xing Yao. "Optimized explicit finite-difference schemes for spatial derivatives using maximum norm." Journal of Computational Physics 250 (2013)

       Conventional fd coefficients formula from:
       3. Chunlei Chu, Paul Stoffa. "Determination of finite-difference weights using scaled binomial windows" Geophysics 77.3 (2012)
       */
    {
      int idim, ii;
#define ndim 3
      float d2[ndim];
      float d1[ndim];
      float *d2_fdcoef;
      float *ccc;
      float *d1_fdcoef;
      float *bbb;

      d2[0] = dz*dz; d2[1] = dx*dx; d2[2] = dy*dy;
      d1[0] = dz; d1[1] = dx; d1[2] = dy;

      if (d_order == 2) {
        d2_fdcoef = sf_floatalloc(ndim*nop+1);
        if (is_optimized) ccc= optimal_fdcoef(nop,2);
        else ccc = normal_fdcoef(nop,d_order);

        d2_fdcoef[0] = 0.f;
        for (idim=0; idim<ndim; idim++) {
          for (ii=1; ii<=nop; ii++) {
            d2_fdcoef[idim*nop+ii] = ccc[ii]/d2[idim];
            d2_fdcoef[0] += d2_fdcoef[idim*nop+ii];
          }
        }
        d2_fdcoef[0] *= - 2.0f;
        return d2_fdcoef;
      } else {
        d1_fdcoef = sf_floatalloc(ndim*nop+1);

        if (is_optimized && nop<=6) bbb = optimal_fdcoef(nop,1);
        else bbb = normal_fdcoef(nop,d_order);
        d1_fdcoef[0] = 0.0f;
        for (idim=0; idim<ndim; idim++) {
          for (ii=1; ii<=nop; ii++) {
            d1_fdcoef[idim*nop+ii] = bbb[ii]/d1[idim];
          }
        }
        return d1_fdcoef;
      }
    }

  int
    factorial(int n)
    {
      int i;
      int result = 1;
      for (i=1; i<=n; ++i)
        result *= i;
      return result;
    }

  float*
    normal_fdcoef(int nop, const int d_order)
    {
      int n;
      float *cc = calloc(nop+1,sizeof(float));
      float *bb = calloc(nop+1,sizeof(float));
      int halfN_fact = factorial(nop);
      halfN_fact *= halfN_fact;
      for (n=1; n<=nop; n++) {
        cc[n] = - 2.f / (n*n) * cos(n*SF_PI) * halfN_fact/factorial(nop+n)/factorial(nop-n);
        bb[n] = cc[n]*n/2.f;
      }
      if (d_order == 1) return bb;
      else return cc;
    }

  float*
    optimal_fdcoef(int nop, const int d_order)
    {
      float *cc;
      float *bb;

      if (d_order == 2)  {
        float *opt_c[11];
        float opt_c1[2] = {0.f, 1.f};
        float opt_c2[3] = {0.f, 1.369074f, - 0.09266816};
        float opt_c3[4] = {0.f, 1.573661f, - 0.1820268f, 0.01728053f};
        float opt_c4[5] = {0.f, 1.700010f, - 0.2554615f, 0.04445392f, - 0.004946851f};
        float opt_c5[6] = {0.f, 1.782836f, - 0.3124513f, 0.07379487f, - 0.01532122f,
          0.001954439f};
        float opt_c6[7] = {0.f, 1.837023f, - 0.3538895f, 0.09978343f, - 0.02815486f,
          0.006556587f, - 0.0009405699f};
        float opt_c7[8] = {0.f, 1.874503f, - 0.3845794f, 0.1215162f, - 0.04121749f,
          0.01295522f, - 0.003313813f, 0.0005310053f};
        float opt_c8[9] = {0.f, 1.901160f, - 0.4074304f, 0.1390909f, - 0.05318775f,
          0.02004823f, - 0.006828249f, 0.001895771f, - 0.0003369052f};
        float opt_c9[10] = {0.f, 1.919909f, - 0.4240446f, 0.1526043f, - 0.06322328f,
          0.02676005f, - 0.01080739f, 0.003907747f, - 0.001158024f,
          0.0002240247f};
        float opt_c10[11] = {0.f, 1.918204f, - 0.4225858f, 0.1514992f, - 0.06249474f,
          0.02637196f, - 0.01066631f, 0.003915625f, - 0.001219872f,
          0.0002863976f, - 0.00003744830f};
        opt_c[1] = opt_c1;  opt_c[2] = opt_c2;  opt_c[3] = opt_c3;  opt_c[4] = opt_c4;
        opt_c[5] = opt_c5;  opt_c[6] = opt_c6;  opt_c[7] = opt_c7;  opt_c[8] = opt_c8;
        opt_c[9] = opt_c9;  opt_c[10] = opt_c10;
        cc = sf_floatalloc(nop+1);
        memcpy(cc,opt_c[nop],sizeof(float)*(nop+1));
        return cc;
      } else {
        float *opt_b[7];
        float opt_b1[2] = {0.0f, 0.5f};
        float opt_b2[3] = {0.0f, 0.67880327, - 0.08962729};
        float opt_b3[4] = {0.0f, 0.77793115, - 0.17388691, 0.02338713};
        float opt_b4[5] = {0.0f, 0.84149635, - 0.24532989, 0.06081891, - 0.00839807};
        float opt_b5[6] = {0.0f, 0.88414717, - 0.30233648, 0.10275057, - 0.02681517,
          0.00398089};
        float opt_b6[7] = {0.0f, 0.91067892, - 0.34187892, 0.13833962, - 0.04880710,
          0.01302148, - 0.00199047};

        opt_b[1] = opt_b1;  opt_b[2] = opt_b2;  opt_b[3] = opt_b3;
        opt_b[4] = opt_b4;  opt_b[5] = opt_b5;  opt_b[6] = opt_b6;
        bb = sf_floatalloc(nop+1);
        memcpy(bb,opt_b[nop],sizeof(float)*(nop+1));
        return bb;
      }
    }

  void
    apply_abc(float*** restrict uu2, float*** restrict uu1, int nz, int nx, int ny, int nbd,
        abcone3d abc, int nop, float* restrict damp)
    {
      int ix, iy, iz, ib;
      float damp_ib;
      float uu2_bc;
      int nxpad = nx + 2*nbd;
      int nypad = ny + 2*nbd;
      int nzpad = nz + 2*nbd;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(ix,iy,iz,ib,damp_ib,uu2_bc)
#endif
      for (iy=0; iy<nypad; iy++) {
        for (ix=0; ix<nxpad; ix++) {
          for (ib=nbd-nop; ib<nbd; ib++) {
            iz = nbd-ib-1;
            uu2[iy][ix][iz] = uu1[iy][ix][iz+1]
              + (uu1[iy][ix][iz] - uu2[iy][ix][iz+1])*abc->bzl[iy][ix];
            iz = nzpad-nbd+ib;
            uu2[iy][ix][iz] = uu1[iy][ix][iz-1]
              + (uu1[iy][ix][iz] - uu2[iy][ix][iz-1])*abc->bzh[iy][ix];
          }
          if (damp != NULL) {
            for (ib=0; ib<nbd-nop; ib++) {
              damp_ib = damp[ib];
              iz = nbd-ib-1;
              uu2_bc = uu1[iy][ix][iz+1] + (uu1[iy][ix][iz] - uu2[iy][ix][iz+1])*abc->bzl[iy][ix];
              uu2[iy][ix][iz] = uu2_bc*(1.f - damp_ib) + uu2[iy][ix][iz]*damp_ib;
              iz = nzpad-nbd+ib;
              uu2_bc = uu1[iy][ix][iz-1] + (uu1[iy][ix][iz] - uu2[iy][ix][iz-1])*abc->bzh[iy][ix];
              uu2[iy][ix][iz] = uu2_bc*(1.f - damp_ib) + uu2[iy][ix][iz]*damp_ib;
            }
          }
          if (abc->free)  memset(uu2[iy][ix], 0, sizeof(float)*(nbd+1));
        }
      }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(ix,iy,iz,ib,damp_ib,uu2_bc)
#endif
      for (iy=0; iy<nypad; iy++) {
        for (iz=0;iz<nzpad; iz++) {
          for (ib=nbd-nop; ib<nbd; ib++) {
            ix = nbd-ib-1;
            uu2[iy][ix][iz] = uu1[iy][ix+1][iz]
              + (uu1[iy][ix][iz] - uu2[iy][ix+1][iz])*abc->bxl[iy][iz];
            ix = nxpad-nbd+ib;
            uu2[iy][ix][iz] = uu1[iy][ix-1][iz]
              + (uu1[iy][ix][iz] - uu2[iy][ix-1][iz])*abc->bxh[iy][iz];
          }
          if (damp != NULL) {
            for (ib=0; ib<nbd-nop; ib++) {
              damp_ib = damp[ib];
              ix = nbd-ib-1;
              uu2_bc = uu1[iy][ix+1][iz] + (uu1[iy][ix][iz] - uu2[iy][ix+1][iz])*abc->bxl[iy][iz];
              uu2[iy][ix][iz] = uu2_bc*(1.f - damp_ib) + uu2[iy][ix][iz]*damp_ib;
              ix = nxpad-nbd+ib;
              uu2_bc = uu1[iy][ix-1][iz] + (uu1[iy][ix][iz] - uu2[iy][ix-1][iz])*abc->bxh[iy][iz];
              uu2[iy][ix][iz] = uu2_bc*(1.f - damp_ib) + uu2[iy][ix][iz]*damp_ib;
            }
          }
        }
      }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(ix,iy,iz,ib,damp_ib,uu2_bc)
#endif
      for (ix=0; ix<nxpad; ix++) {
        for (iz=0;iz<nzpad; iz++) {
          for (ib=nbd-nop; ib<nbd; ib++) {
            iy = nbd-ib-1;
            uu2[iy][ix][iz] = uu1[iy+1][ix][iz]
              + (uu1[iy][ix][iz] - uu2[iy+1][ix][iz])*abc->byl[ix][iz];
            iy = nypad-nbd+ib;
            uu2[iy][ix][iz] = uu1[iy-1][ix][iz]
              + (uu1[iy][ix][iz] - uu2[iy-1][ix][iz])*abc->byh[ix][iz];
          }
          if (damp != NULL) {
            for (ib=0; ib<nbd-nop; ib++) {
              damp_ib = damp[ib];
              iy = nbd-ib-1;
              uu2_bc = uu1[iy+1][ix][iz] + (uu1[iy][ix][iz] - uu2[iy+1][ix][iz])*abc->byl[ix][iz];
              uu2[iy][ix][iz] = uu2_bc*(1.f - damp_ib) + uu2[iy][ix][iz]*damp_ib;
              iy = nypad-nbd+ib;
              uu2_bc = uu1[iy-1][ix][iz] + (uu1[iy][ix][iz] - uu2[iy-1][ix][iz])*abc->byh[ix][iz];
              uu2[iy][ix][iz] = uu2_bc*(1.f - damp_ib) + uu2[iy][ix][iz]*damp_ib;
            }
          }
        }
      }
      return;
    }

