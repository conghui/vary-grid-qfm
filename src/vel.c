#include "vel.h"

#include <rsf.h>

static float g_w0;
static float g_qfact;
static float g_gamma;

vel_t *new_vel(float ***vel, int nz, int nx, int ny,
    float oz, float ox, float oy,
    float dz, float dx, float dy,
    float w0, float qfact) 
{
  vel_t *v = malloc(sizeof *v);
  v->nz = nz; v->nx = nx; v->ny = ny;
  v->oz = oz; v->ox = ox; v->oy = oy;
  v->dz = dz; v->dx = dx; v->dy = dy;
  v->dat = vel;
  v->vgamma = NULL;

  g_w0 = w0; g_qfact= qfact; 
  g_gamma = 1.0/SF_PI*atan(2*SF_PI/qfact);
  return v;
}
