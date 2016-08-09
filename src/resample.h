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

#endif /* end of include guard: RESAMPLE_H_WQCWR39E */
