#include "vel.h"

#include <rsf.h>

static float gs_w0;
static float gs_qfact;
static float gs_gamma;

vel_t *clone_vel(float ***vel, int nz, int nx, int ny,
    float oz, float ox, float oy,
    float dz, float dx, float dy,
    float w0, float qfact)
{
  vel_t *v = malloc(sizeof *v);
  v->n1 = nz; v->n2 = nx; v->n3 = ny;
  v->o1 = oz; v->o2 = ox; v->o3 = oy;
  v->d1 = dz; v->d2 = dx; v->d3 = dy;
  v->dat = vel;
  v->vgamma = NULL;

  gs_w0 = w0; gs_qfact= qfact;
  gs_gamma = 1.0/SF_PI*atan(2*SF_PI/qfact);
  return v;
}


void vmin_vmax_dmin_dmax(const vel_t *vel, float *vmin, float *vmax, float *dmin, float *dmax)
{
  *vmin = *dmin = 99999999;
  *vmax = *dmax = 0;

  for (int i3 = 0; i3 < vel->n3; i3++) {
    for (int i2 = 0; i2 < vel->n2; i2++) {
      for (int i1 = 0; i1 < vel->n1; i1++) {
        *vmin = fminf(*vmin, vel->dat[i3][i2][i1]);
        *vmax = fmaxf(*vmax, vel->dat[i3][i2][i1]);
      }
    }
  }

  *dmin = fminf(fminf(vel->d1, vel->d2), vel->d3);
  *dmax = fmaxf(fmaxf(vel->d1, vel->d2), vel->d3);
}
