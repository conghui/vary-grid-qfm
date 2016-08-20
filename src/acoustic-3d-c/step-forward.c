#include "step-forward.h"

#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
#endif

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
void step_forward(float*** restrict u0, float*** restrict u1,
    float*** restrict vel, float*** restrict rho,
    float* restrict fdcoef_d2, float* restrict fdcoef_d1,
    int nop, int nzpad, int nxpad, int nypad)
{
  int ix, iy, iz, iop;
  float c0 = fdcoef_d2[0];
  float *cz = &fdcoef_d2[0];
  float *cx = &fdcoef_d2[nop];
  float *cy = &fdcoef_d2[nop+nop];
  float *bz, *bx, *by;
  float drho_dot_du;
  float du_z = 0.f, du_x = 0.f, du_y = 0.f;
  float drho_z = 0.f, drho_x = 0.f, drho_y = 0.f;
  float lap;

  if (rho != NULL) { /* variable density */
    bz = &fdcoef_d1[0];
    bx = &fdcoef_d1[nop];
    by = &fdcoef_d1[nop+nop];
  }

#ifdef __SSE__
#define vType __m128
#define _prefix_loadu_ps _mm_loadu_ps
#define _prefix_storeu_ps _mm_storeu_ps
#define _prefix_add_ps _mm_add_ps
#define _prefix_sub_ps _mm_sub_ps
#define _prefix_mul_ps _mm_mul_ps
#define _prefix_div_ps _mm_div_ps
#define _prefix_set1_ps _mm_set1_ps
#elif defined __AVX__
#define vType __m256
#define _prefix_loadu_ps _mm256_loadu_ps
#define _prefix_storeu_ps _mm256_storeu_ps
#define _prefix_add_ps _mm256_add_ps
#define _prefix_sub_ps _mm256_sub_ps
#define _prefix_mul_ps _mm256_mul_ps
#define _prefix_div_ps _mm256_div_ps
#define _prefix_set1_ps _mm256_set1_ps
#endif

#ifdef _OPENMP
#pragma omp parallel for				\
  schedule(static,1)					\
  shared(nxpad,nypad,nzpad,u0,u1,vel,c0,cx,cy,cz) \
  private(ix,iy,iz,iop,drho_dot_du,du_z,du_x,du_y,drho_z,drho_x,drho_y,lap)
#endif
#if defined (__SSE__) || defined (__AVX__)
  for (iy=nop; iy<nypad-nop; iy++) {
    for (ix=nop; ix<nxpad-nop; ix++) {
#ifdef __SSE__
      for (iz=nop; iz<nzpad-nop; iz+=4) {
#elif defined __AVX__
        for (iz=nop; iz<nzpad-nop; iz+=8) {
#endif
          /* load u0 u1 vel from memory to register */
          vType vec_u0 = _prefix_loadu_ps(&u0[iy][ix][iz]);
          vType vec_u1 = _prefix_loadu_ps(&u1[iy][ix][iz]);
          vec_u0 = _prefix_sub_ps(_prefix_mul_ps(_prefix_set1_ps(2.f),vec_u1),vec_u0);

          vType vec_lap = _prefix_mul_ps(vec_u1,_prefix_set1_ps(c0));
          for (iop=1; iop<=nop; iop++) {
            /* z axis derivative <1st axis> */
            vec_lap = _prefix_add_ps(vec_lap,
                _prefix_mul_ps(
                  _prefix_set1_ps(cz[iop]),
                  _prefix_add_ps(
                    _prefix_loadu_ps(&u1[iy][ix][iz+iop]),
                    _prefix_loadu_ps(&u1[iy][ix][iz-iop])
                    )
                  )
                );
            /* x axis derivative <2nd axis> */
            vec_lap = _prefix_add_ps(vec_lap,
                _prefix_mul_ps(
                  _prefix_set1_ps(cx[iop]),
                  _prefix_add_ps(
                    _prefix_loadu_ps(&u1[iy][ix+iop][iz]),
                    _prefix_loadu_ps(&u1[iy][ix-iop][iz])
                    )
                  )
                );
            /* y axis derivative <3rd axis> */
            vec_lap = _prefix_add_ps(vec_lap,
                _prefix_mul_ps(
                  _prefix_set1_ps(cy[iop]),
                  _prefix_add_ps(
                    _prefix_loadu_ps(&u1[iy+iop][ix][iz]),
                    _prefix_loadu_ps(&u1[iy-iop][ix][iz])
                    )
                  )
                );
          }

          if (rho != NULL) {
            vType vec_du_z = _prefix_set1_ps(0.f);
            vType vec_du_x = _prefix_set1_ps(0.f);
            vType vec_du_y = _prefix_set1_ps(0.f);
            vType vec_drho_z = _prefix_set1_ps(0.f);
            vType vec_drho_x = _prefix_set1_ps(0.f);
            vType vec_drho_y = _prefix_set1_ps(0.f);
            for (iop=1; iop<=nop; iop++) {
              vec_du_z = _prefix_add_ps(vec_du_z,
                  _prefix_mul_ps(
                    _prefix_set1_ps(bz[iop]),
                    _prefix_sub_ps(
                      _prefix_loadu_ps(&u1[iy][ix][iz+iop]),
                      _prefix_loadu_ps(&u1[iy][ix][iz-iop])
                      )
                    )
                  );
              vec_drho_z = _prefix_add_ps(vec_drho_z,
                  _prefix_mul_ps(
                    _prefix_set1_ps(bz[iop]),
                    _prefix_sub_ps(
                      _prefix_loadu_ps(&rho[iy][ix][iz+iop]),
                      _prefix_loadu_ps(&rho[iy][ix][iz-iop])
                      )
                    )
                  );
              vec_du_x = _prefix_add_ps(vec_du_x,
                  _prefix_mul_ps(
                    _prefix_set1_ps(bx[iop]),
                    _prefix_sub_ps(
                      _prefix_loadu_ps(&u1[iy][ix+iop][iz]),
                      _prefix_loadu_ps(&u1[iy][ix-iop][iz])
                      )
                    )
                  );
              vec_drho_x = _prefix_add_ps(vec_drho_x,
                  _prefix_mul_ps(
                    _prefix_set1_ps(bx[iop]),
                    _prefix_sub_ps(
                      _prefix_loadu_ps(&rho[iy][ix+iop][iz]),
                      _prefix_loadu_ps(&rho[iy][ix-iop][iz])
                      )
                    )
                  );
              vec_du_y = _prefix_add_ps(vec_du_y,
                  _prefix_mul_ps(
                    _prefix_set1_ps(by[iop]),
                    _prefix_sub_ps(
                      _prefix_loadu_ps(&u1[iy+iop][ix][iz]),
                      _prefix_loadu_ps(&u1[iy-iop][ix][iz])
                      )
                    )
                  );
              vec_drho_y = _prefix_add_ps(vec_drho_y,
                  _prefix_mul_ps(
                    _prefix_set1_ps(by[iop]),
                    _prefix_sub_ps(
                      _prefix_loadu_ps(&rho[iy+iop][ix][iz]),
                      _prefix_loadu_ps(&rho[iy-iop][ix][iz])
                      )
                    )
                  );
            }
            vec_lap = _prefix_sub_ps(vec_lap,
                _prefix_div_ps(
                  _prefix_add_ps(
                    _prefix_mul_ps(vec_du_z,vec_drho_z),
                    _prefix_add_ps(
                      _prefix_mul_ps(vec_du_x,vec_drho_x),
                      _prefix_mul_ps(vec_du_y,vec_drho_y)
                      )
                    ),
                  _prefix_loadu_ps(&rho[iy][ix][iz])
                  )
                );
          }
          vec_u0 = _prefix_add_ps(vec_u0,_prefix_mul_ps(_prefix_loadu_ps(&vel[iy][ix][iz]),vec_lap));
          _prefix_storeu_ps(&u0[iy][ix][iz],vec_u0);
        }
      }
    }
#undef vType
#undef _prefix_loadu_ps
#undef _prefix_storeu_ps
#undef _prefix_add_ps
#undef _prefix_sub_ps
#undef _prefix_mul_ps
#undef _prefix_div_ps
#undef _prefix_set1_ps

#else
    for (iy=nop; iy<nypad-nop; iy++) {
      for (ix=nop; ix<nxpad-nop; ix++) {
        for (iz=nop; iz<nzpad-nop; iz++) {
          lap = u1[iy][ix][iz]*c0;
          for (iop=1; iop<=nop; iop++) {
            lap +=  (u1[iy][ix][iz-iop] + u1[iy][ix][iz+iop]) * cz[iop]
                  + (u1[iy][ix-iop][iz] + u1[iy][ix+iop][iz]) * cx[iop]
                  + (u1[iy-iop][ix][iz] + u1[iy+iop][ix][iz]) * cy[iop];
          }
          if (rho != NULL) { /* variable density term */
            du_z = du_x = du_y = drho_z = drho_x = drho_y = 0.f;
            for (iop=1; iop<=nop; iop++) {
              du_z += (u1[iy][ix][iz+iop] - u1[iy][ix][iz-iop]) * bz[iop];
              du_x += (u1[iy][ix+iop][iz] - u1[iy][ix-iop][iz]) * bx[iop];
              du_y += (u1[iy+iop][ix][iz] - u1[iy-iop][ix][iz]) * by[iop];
              drho_z += (rho[iy][ix][iz+iop] - rho[iy][ix][iz-iop]) * bz[iop];
              drho_x += (rho[iy][ix+iop][iz] - rho[iy][ix-iop][iz]) * bx[iop];
              drho_y += (rho[iy+iop][ix][iz] - rho[iy-iop][ix][iz]) * by[iop];
            }
            drho_dot_du = (du_z*drho_z + du_x*drho_x + du_y*drho_y)/rho[iy][ix][iz];
            lap -= drho_dot_du;
          }
          u0[iy][ix][iz] = 2.*u1[iy][ix][iz] - u0[iy][ix][iz] + vel[iy][ix][iz]*lap;
        }
      }
    }
#endif
    return;
}

void step_forward_q(float*** restrict u0, float*** restrict u1,
    float*** restrict vel, float *** restrict vgamma, float*** restrict rho,
    float* restrict fdcoef_d2, float* restrict fdcoef_d1,
    int nop, int nzpad, int nxpad, int nypad)
{
  float c0 = fdcoef_d2[0];
  float *cz = &fdcoef_d2[0];
  float *cx = &fdcoef_d2[nop];
  float *cy = &fdcoef_d2[nop+nop];

  (void)(rho);
  (void)(fdcoef_d1);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int iy=nop; iy<nypad-nop; iy++) {
    for (int ix=nop; ix<nxpad-nop; ix++) {
      for (int iz=nop; iz<nzpad-nop; iz++) {
        float tau = vgamma[iy][ix][iz];
        float taucur = 1 - tau;
        float tauprev = tau;
        float lap = (taucur*u1[iy][ix][iz] + tauprev*u0[iy][ix][iz])*c0;

        for (int iop=1; iop<=nop; iop++) {
          lap +=
              (taucur*(u1[iy][ix][iz-iop] + u1[iy][ix][iz+iop]) + tauprev*(u0[iy][ix][iz-iop] + u0[iy][ix][iz+iop])) * cz[iop]
            + (taucur*(u1[iy][ix-iop][iz] + u1[iy][ix+iop][iz]) + tauprev*(u0[iy][ix-iop][iz] + u0[iy][ix+iop][iz])) * cx[iop]
            + (taucur*(u1[iy-iop][ix][iz] + u1[iy+iop][ix][iz]) + tauprev*(u0[iy-iop][ix][iz] + u0[iy+iop][ix][iz])) * cy[iop];
        }
        u0[iy][ix][iz] = 2.*u1[iy][ix][iz] - u0[iy][ix][iz] + vel[iy][ix][iz]*lap;
      }
    }
  }
}
