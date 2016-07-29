#ifndef BOX_H_MLWKBDTF
#define BOX_H_MLWKBDTF

typedef struct {
  int n1, n2, n3, n4, n5;
  float d1, d2, d3, d4, d5;
  float o1, o2, o3, o4, o5;

  float *****val;
} times_t;

times_t *read_times();

#endif /* end of include guard: BOX_H_MLWKBDTF */
