#ifndef VEL_H_BLY3XMED
#define VEL_H_BLY3XMED

#include "model.h"

typedef struct {
  int   n1, n2, n3;
  float o1, o2, o3;
  float d1, d2, d3;
  float ***dat;
  float ***vgamma;
} vel_t;

void vmin_vmax_dmin_dmax(const vel_t *vel, float *vmin, float *vmax, float *dmin, float *dmax);

vel_t *clone_vel(float ***vel, int nz, int nx, int ny,
    float oz, float ox, float oy,
    float dz, float dx, float dy,
    float w0, float qfact);

vel_t *read_vel(const char *tag);
modeling_t make_modeling(const vel_t *v);

#endif /* end of include guard: VEL_H_BLY3XMED */
