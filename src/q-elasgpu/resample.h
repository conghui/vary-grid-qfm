#ifndef RESAMPLE_H_WQCWR39E
#define RESAMPLE_H_WQCWR39E

#include <rsf.h>

void init_sinc_table(int nsinc, int npts);
void interpfield_(float ***oldf, float ***newf, bool extend,
    int on1, float oo1, float od1,  /* old */
    int on2, float oo2, float od2,
    int on3, float oo3, float od3,
    int nn1, float no1, float nd1,  /* new */
    int nn2, float no2, float nd2,
    int nn3, float no3, float nd3);

void interp_den_vel_(
    float ***full_h_ro, float ***full_h_c11, float ***full_h_c22, float ***full_h_c33, 
    float ***full_h_c44, float ***full_h_c55, float ***full_h_c66, float ***full_h_c12, 
    float ***full_h_c13, float ***full_h_c23, 
    float ***h_ro, float ***h_c11, float ***h_c22, float ***h_c33, float ***h_c44, 
    float ***h_c55, float ***h_c66, float ***h_c12, float ***h_c13, float ***h_c23,
    int on1, float oo1, float od1,  /* old */
    int on2, float oo2, float od2,
    int on3, float oo3, float od3,
    int nn1, float no1, float nd1,  /* new */
    int nn2, float no2, float nd2,
    int nn3, float no3, float nd3);

void interp_wavefield_(
    float ***o_umx, float ***o_uox,  float ***o_umy,  float ***o_uoy,  float ***o_umz,  float ***o_uoz,
    float ***n_umx, float ***n_uox,  float ***n_umy,  float ***n_uoy,  float ***n_umz,  float ***n_uoz,
    int on1, float oo1, float od1,  /* old */
    int on2, float oo2, float od2,
    int on3, float oo3, float od3,
    int nn1, float no1, float nd1,  /* new */
    int nn2, float no2, float nd2,
    int nn3, float no3, float nd3);

#endif /* end of include guard: RESAMPLE_H_WQCWR39E */
