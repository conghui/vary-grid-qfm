#ifndef RESAMPLE_H_WQCWR39E
#define RESAMPLE_H_WQCWR39E

#include "model.h"
#include <rsf.h>

void init_sinc_table(int nsinc, int npts);
void interpfield(const modeling_t *olds, const modeling_t *news, float ***oldf, float ***newf, bool extend);

#endif /* end of include guard: RESAMPLE_H_WQCWR39E */
