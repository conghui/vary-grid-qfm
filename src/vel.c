#include "vel.h"
#include "resample.h"

#include <rsf.h>

static float gs_w0;
static float gs_qfact;
static float gs_gamma;

static void interpfield(const modeling_t *olds, const modeling_t *news, float ***oldf, float ***newf, bool extend)
{
  interpfield_(oldf, newf, extend,
    olds->n1, olds->o1, olds->d1,  /* old */
    olds->n2, olds->o2, olds->d2,
    olds->n3, olds->o3, olds->d3,
    news->n1, news->o1, news->d1,  /* new */
    news->n2, news->o2, news->d2,
    news->n3, news->o3, news->d3);
}

modeling_t make_modeling(const vel_t *v)
{
  modeling_t m;
  m.n1 = v->n1; m.o1 = v->o1; m.d1 = v->d1;
  m.n2 = v->n2; m.o2 = v->o2; m.d2 = v->d2;
  m.n3 = v->n3; m.o3 = v->o3; m.d3 = v->d3;

  return m;
}

vel_t *read_vel(const char *tag)
{
  sf_file fvel = sf_input(tag);
  sf_seek(fvel, 0, SEEK_SET);

  vel_t *v = malloc(sizeof *v);
  sf_histint(fvel, "n1", &v->n1);
  sf_histint(fvel, "n2", &v->n2);
  sf_histint(fvel, "n3", &v->n3);
  sf_histfloat(fvel, "o1", &v->o1);
  sf_histfloat(fvel, "o2", &v->o2);
  sf_histfloat(fvel, "o3", &v->o3);
  sf_histfloat(fvel, "d1", &v->d1);
  sf_histfloat(fvel, "d2", &v->d2);
  sf_histfloat(fvel, "d3", &v->d3);

  v->dat = sf_floatalloc3(v->n1, v->n2, v->n3);
  sf_floatread(v->dat[0][0], v->n1 * v->n2 * v->n3, fvel);

  return v;
}

vel_t *clone_vel(float ***vel, int nz, int nx, int ny,
    float oz, float ox, float oy,
    float dz, float dx, float dy,
    float w0, float qfact)
{
  vel_t *v = malloc(sizeof *v);
  v->n1 = nz; v->n2 = nx; v->n3 = ny;
  v->o1 = oz; v->o2 = ox; v->o3 = oy;
  v->d1 = dz; v->d2 = dx; v->d3 = dy;
  v->dat = sf_floatalloc3(nz, nx, ny);
  memcpy(v->dat[0][0], vel[0][0], nz*nx*ny*sizeof(float));
  v->vgamma = sf_floatalloc3(1,1,1);

  gs_w0 = w0; gs_qfact= qfact;
  gs_gamma = 1.0/SF_PI*atan(2*SF_PI/qfact);

  sf_warning("gs_gamma: %f", gs_gamma);
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

void resample_vel(const modeling_t *olds, const modeling_t *news, const vel_t *oldv, vel_t *newv) {
  free(**newv->dat); free(*newv->dat); free(newv->dat);
  free(**newv->vgamma); free(*newv->vgamma); free(newv->vgamma);

  newv->dat = sf_floatalloc3(news->n1, news->n2, news->n3);
  newv->vgamma = sf_floatalloc3(news->n1, news->n2, news->n3);

  interpfield(olds, news, oldv->dat, newv->dat, true);

  for (int i3 = 0; i3 < news->n3; i3++) {
    for (int i2 = 0; i2 < news->n2; i2++) {
      for (int i1 = 0; i1 < news->n1; i1++) {
        newv->vgamma[i3][i2][i1] = -powf(newv->dat[i3][i2][i1],2*gs_gamma-1) * powf(gs_w0, 2 * gs_gamma) * sin(SF_PI * gs_gamma) / news->dt;
      }
    }
  }
}

void resample_p(const modeling_t *olds, const modeling_t *news, float ****p)
{
  float ***newp = sf_floatalloc3(news->n1, news->n2, news->n3);
  interpfield(olds, news, *p, newp, false);

  free((*p)[0][0]); free((*p)[0]); free(*p);
  *p = newp;
}
