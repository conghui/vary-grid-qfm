#include <rsf.h>
#include "box.h"

static int   gs_timeblocks;
static float gs_vmin;
static float gs_vmax;
static float gs_dmin;
static float gs_dmax;
static float gs_dtstable; // maximum sampling in time for stability
static float gs_dtbig;
static float gs_maxf; // maximum frequency
static int   gs_nb; // boundary
static int   gs_error;
static float gs_errorfact;
static float gs_qfact;
static float gs_downfact;

static float calgoodsampling(dt) {
  return 0.95 * dt;
}

static void find_exts(const int *ns, const float *os, const float *ds, float ***minT, float maxT, float *mymin, float *mymax)
{
  const int D = 3; // # of dimensions
  for (int i = 0; i < D; i++) {
    mymin[i] = os[i] + ns[i] * ds[i];
    mymax[i] = os[i];
  }

  for (int i3 = 0; i3 < ns[2]; i3++) {
    for (int i2 = 0; i2 < ns[1]; i2++) {
      for (int i1 = 0; i1 < ns[0]; i1++) {
        if (minT[i3][i2][i1] <= maxT) {
          mymin[0] = fminf(mymin[0], os[0] + ds[0] * i1);
          mymin[1] = fminf(mymin[1], os[1] + ds[1] * i2);
          mymin[2] = fminf(mymin[2], os[2] + ds[2] * i3);
          mymax[0] = fmaxf(mymax[0], os[0] + ds[0] * i1);
          mymax[1] = fmaxf(mymax[1], os[1] + ds[1] * i2);
          mymax[2] = fmaxf(mymax[2], os[2] + ds[2] * i3);
        }
      }
    }
  }
}

static void create_modeling(modeling_t *mod, const times_t *times, float ***minT, float vmin, float dmax, float maxf, float timemin, float timemax, float qfact, float downfact, float errorfact)
{
  int   ns[3] = {times->n1, times->n2, times->n3};
  float os[3] = {times->o1, times->o2, times->o3};
  float ds[3] = {times->d1, times->d2, times->d3};
  float mymin[3], mymax[3];
  float dsamp;

  for (int i = 0; i < 3; i++) {
    sf_warning("ns[%d]: %d, os[%d]: %f, ds[%d]: %f", i, ns[i], i, os[i], i, ds[i]);
  }
  sf_warning("timemax: %f", timemax);

  find_exts(ns, os, ds, minT, timemax, mymin, mymax);

  for (int i = 0; i < 3; i++) {
    sf_warning("mymin[%d]: %f, mymax[%d]: %f", i, mymin[i], i, mymax[i]);
  }

  if (timemin > 0.0001) {
    sf_warning("compute __dsamp");
    float ff = 1.0 - 2.0 * SF_PI / qfact;
    float cycleskill = log(downfact) / log(ff);
    float fkill = cycleskill / fmaxf(timemin, 0.001);
    dsamp = fmax(dmax, vmin/fminf(maxf, fkill) / 3.3 / errorfact);
  } else {
    sf_warning("direct set __dsamp to dmax");
    dsamp = dmax;
  }

  sf_warning("dsamp: %f", dsamp);
  /*exit(0);*/
}

times_t *read_times() {
  times_t *t = malloc(sizeof *t);

  sf_file file_times = sf_input("times");
  sf_histint(file_times, "n1", &t->n1);
  sf_histint(file_times, "n2", &t->n2);
  sf_histint(file_times, "n3", &t->n3);
  sf_histint(file_times, "n4", &t->n4);
  sf_histint(file_times, "n5", &t->n5);
  sf_histfloat(file_times, "o1", &t->o1);
  sf_histfloat(file_times, "o2", &t->o2);
  sf_histfloat(file_times, "o3", &t->o3);
  sf_histfloat(file_times, "o4", &t->o4);
  sf_histfloat(file_times, "o5", &t->o5);
  sf_histfloat(file_times, "d1", &t->d1);
  sf_histfloat(file_times, "d2", &t->d2);
  sf_histfloat(file_times, "d3", &t->d3);
  sf_histfloat(file_times, "d4", &t->d4);
  sf_histfloat(file_times, "d5", &t->d5);

  t->val = sf_floatalloc5(t->n1, t->n2, t->n3, t->n4, t->n5);

  int size = t->n1 * t->n2 * t->n3 * t->n4 * t->n5;
  sf_floatread(t->val[0][0][0][0], size, file_times);

  return t;
}

void init_box(int   timeblocks, float vmin, float vmax, float dmin, float dmax, float maxf, int   nb, int   error, float errorfact, float qfact, float downfact)
{
  gs_timeblocks = timeblocks;
  gs_vmin       = vmin;
  gs_vmax       = vmax;
  gs_dmin       = dmin;
  gs_dmax       = dmax;
  gs_maxf       = maxf;
  gs_nb         = nb;
  gs_error      = error;
  gs_errorfact  = errorfact;
  gs_qfact      = qfact;
  gs_downfact   = downfact;

  gs_dtstable   = 0.499 * dmin / vmax;
  gs_dtbig      = calgoodsampling(gs_dtstable);
}

void calc_shot_box(const times_t *times, const pt3d *src3d, const pt3d *rec3d, int nr, int nt, float dt) {

  int i4 = (src3d->x - times->o4) / times->d4 + 0.5;
  int i5 = (src3d->y - times->o5) / times->d5 + 0.5;

  /*sf_warning("i4: %d, i5: %d", i4, i5);*/

  float ***minT = sf_floatalloc3(times->n1, times->n2, times->n3);
  memcpy(minT[0][0], times->val[i5][i4][0][0],
      sizeof(float) * times->n1 * times->n2 * times->n3);

  for (int i = 0; i < nr; i++) {
    int i5 = fminf(fmaxf(round((rec3d[i].y - times->o5) / times->d5 + 0.5), 0), times->n5);
    int i4 = fminf(fmaxf(round((rec3d[i].x - times->o4) / times->d4 + 0.5), 0), times->n4);

    /*sf_warning("i4, i5: %d, %d", i4, i5);*/
    for (int i3 = 0; i3 < times->n3; i3++) {
      for (int i2 = 0; i2 < times->n2; i2++) {
        for (int i1 = 0; i1 < times->n1; i1++) {
          minT[i3][i2][i1] = fminf(minT[i3][i2][i1], times->val[i5][i4][i3][i2][i1]);
        }
      }
    }
  }

  sf_file file_minT = sf_output("minT");
  sf_putint(file_minT, "n1", times->n1);
  sf_putint(file_minT, "n2", times->n2);
  sf_putint(file_minT, "n3", times->n3);
  sf_floatwrite(minT[0][0], times->n1 * times->n2 * times->n3, file_minT);

  box_t domain;
  int timeblocks = gs_timeblocks; // global variable
  domain.hyper = malloc(timeblocks * sizeof(*domain.hyper));
  /// time fdm3d doesn't initialized
  float blocktime = nt * dt / timeblocks;

  sf_warning("timeblocks: %d", timeblocks);
  sf_warning("blocktime: %f", blocktime);

  for (int i = 0; i < timeblocks; i++) {
    float timemin = i * blocktime;
    float timemax = timemin + blocktime;

    /*create_modeling(&domain.hyper[i], times, minT, timemax);*/
  create_modeling(&domain.hyper[i], times, minT, gs_vmin, gs_dmax, gs_maxf, timemin, timemax, gs_qfact, gs_downfact, gs_errorfact);

  }
}
