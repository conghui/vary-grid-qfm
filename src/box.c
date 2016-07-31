#include <rsf.h>
#include "box.h"

static int   gs_timeblocks;
static float gs_vmin;
static float gs_vmax;
static float gs_dmin;
static float gs_dmax;
static float gs_dtbig;
static float gs_maxf; // maximum frequency
static int   gs_nb; // boundary
static float gs_error;
static float gs_errorfact;
static float gs_qfact;
static float gs_downfact;

static float calgoodsampling(float dt) {
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

static void create_modeling(modeling_t *mod, const vel_t *vel, const times_t *times, float ***minT, float vmin, float vmax, float dmax, float maxf, float timemin, float timemax, float qfact, float downfact, float errorfact, float error, float *totallargecells, float *totalsmallcells)
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

  sf_warning("dmax: %f", dmax);
  if (timemin > 0.0001) {
    sf_warning("compute __dsamp");
    float ff = 1.0 - 2.0 * SF_PI / qfact;
    float cycleskill = log(downfact) / log(ff);
    float fkill = cycleskill / fmaxf(timemin, 0.001);
    /*dsamp = fmax(dmax, vmin/fminf(maxf, fkill) / 3.3 / errorfact);*/
    /// if your unit is: km/s, then divide by 1000
    dsamp = fmax(dmax, vmin/fminf(maxf, fkill) / 3.3 / errorfact / 1000);
  } else {
    sf_warning("direct set __dsamp to dmax");
    dsamp = dmax;
  }

  sf_warning("dsamp: %f", dsamp);

  sf_warning("vel->o1: %f", vel->o1);
  sf_warning("vel->n1: %d", vel->n1);

  float bv[3] = {vel->o1, vel->o2, vel->o3};
  float ev[3] = {vel->o1 + (vel->n1 - 1) * vel->d1,
                 vel->o2 + (vel->n2 - 1) * vel->d2,
                 vel->o3 + (vel->n3 - 1) * vel->d3};

  for (int i = 0; i < 3; i++) {
    sf_warning("bv[%d]: %f, ev[%d]: %f", i, bv[i], i, ev[i]);
  }

  sf_warning("ev_: %f, %f, %f", ev[0], ev[1], ev[2]);
  sf_warning("mymin_: %f, %f, %f", mymin[0], mymin[1], mymin[2]);
  sf_warning("mymax_: %f, %f, %f", mymax[0], mymax[1], mymax[2]);
  sf_warning("error: %f", error);
  sf_warning("dsamp: %f", dsamp);
  float o[3];
  int   n[3];
  for (int i = 0; i < 3; i++) {
    o[i] = fmaxf(bv[i], mymin[i] - dsamp * error);
    mymax[i] = fminf(ev[i], mymax[i] + dsamp * error);
    n[i] = ceilf((mymax[i] - mymin[i]) / dsamp) + 1;
  }

  sf_warning("o(:): %f, %f, %f", o[0], o[1], o[2]);
  sf_warning("mymax(:): %f, %f, %f", mymax[0], mymax[1], mymax[2]);
  sf_warning("n(:): %d, %d, %d", n[0], n[1], n[2]);

  mod->d1 = mod->d2 = mod->d3 = dsamp;
  mod->o1 = o[0] - mod->d1 * mod->nb;
  mod->o2 = o[1] - mod->d2 * mod->nb;
  mod->o3 = o[2] - mod->d3 * mod->nb;
  mod->n1 = n[0] + 2 * mod->nb;
  mod->n2 = n[1] + 2 * mod->nb;
  mod->n3 = n[2] + 2 * mod->nb;

  float dtmax = 0.49 * dsamp / vmax;
  mod->dt = calgoodsampling(dtmax);
  mod->dt = 0.004; // TODO: update dt
  mod->ntblock = (timemax - timemin) / mod->dt;
  mod->dtextra = (timemax - timemin) - mod->dt * mod->ntblock;

  sf_warning("mod->n(:): %d, %d, %d", mod->n1, mod->n2, mod->n3);
  sf_warning("mod->o(:): %f, %f, %f", mod->o1, mod->o2, mod->o3);
  sf_warning("mod->d(:): %f, %f, %f", mod->d1, mod->d2, mod->d3);
  sf_warning("mod->ntblock: %d", mod->ntblock);
  sf_warning("mod->dtextra: %f", mod->dtextra);

  float smallsize = 0.001 * mod->n1 * mod->n2 * mod->n3;
  float fullsize =  0.001*(vel->n1+2*mod->nb)*(vel->n2+2*mod->nb)*(vel->n3+2*mod->nb);
  float small = smallsize * mod->ntblock;
  float large = fullsize * (timemax - timemin) / gs_dtbig;
  /*sf_warning("smallsize_: %f, ntblock: %d, small: %f", smallsize, mod->ntblock, small);*/
  /*sf_warning("small_: %d", small);*/
  *totalsmallcells += small;
  *totallargecells += large;

  sf_warning("dtbig: %f", gs_dtbig);
  sf_warning("smallsize: %f, largesize: %f", smallsize, fullsize);
  sf_warning("small: %f, large: %f", small, large);
  sf_warning("speedup factor: %f", large / small);

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

void init_box(int   timeblocks, float vmin, float vmax, float dmin, float dmax, float maxf, int   nb, float   error, float errorfact, float qfact, float downfact)
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

  float dstable   = 0.499 * dmin / vmax;
  gs_dtbig      = calgoodsampling(dstable);

  sf_warning("dmin: %f, vmax: %f", dmin, vmax);
  sf_warning("dstable: %f", dstable);
  sf_warning("gs_dtbig: %f", gs_dtbig);
}

box_t *calc_shot_box(const vel_t *vel, const times_t *times, const pt3d *src3d, const pt3d *rec3d, int nr, int nt, float dt) {

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

  box_t *domain = malloc(sizeof *domain);
  int timeblocks = gs_timeblocks; // global variable
  domain->timeblocks = timeblocks;
  domain->hyper = malloc(timeblocks * sizeof(*domain->hyper));
  /// time fdm3d doesn't initialized
  float blocktime = nt * dt / timeblocks;

  sf_warning("timeblocks: %d", timeblocks);
  sf_warning("blocktime: %f", blocktime);

  float totallargecells = 0;
  float totalsmallcells = 0;

  for (int i = 0; i < timeblocks; i++) {
    domain->hyper[i].nb = gs_nb;
    float timemin = i * blocktime;
    float timemax = timemin + blocktime;

  create_modeling(&domain->hyper[i], vel, times, minT, gs_vmin, gs_vmax, gs_dmax, gs_maxf, timemin, timemax, gs_qfact, gs_downfact, gs_errorfact, gs_error, &totallargecells, &totalsmallcells);

  }

  float totalspeedup = totallargecells / totalsmallcells;
  sf_warning("timeblocks: %d, totalspeedup: %f", timeblocks, totalspeedup);

  return domain;
}
