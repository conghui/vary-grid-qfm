#ifndef BOX_H_MLWKBDTF
#define BOX_H_MLWKBDTF

#include <rsf.h>
#include "vel.h"
#include "model.h"

typedef struct {
  int n1, n2, n3, n4, n5;
  float d1, d2, d3, d4, d5;
  float o1, o2, o3, o4, o5;

  float *****val;
} times_t;


typedef struct  {
  int nt;
  int totImaging;
  float dt;
  float blockSize;  // propagation time in each block
  int timeblocks;
  modeling_t *hyper;
} box_t;

times_t *read_times();

box_t *calc_shot_box(const vel_t *vel, const times_t *times, const pt3d *src3d, const pt3d *rec3d, int nr, int nt, float dt);

void init_box(int timeblocks, float vmin, float vmax, float dmin, float dmax, float maxf, int nb, float error, float errorfact, float qfact, float downfact, float dt);
#endif /* end of include guard: BOX_H_MLWKBDTF */
