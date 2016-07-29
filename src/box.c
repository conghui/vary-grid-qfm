#include <rsf.h>
#include "box.h"

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
  sf_floatread(t->val[0][0][0][0], size * sizeof(float), file_times);

  return t;
}

/*void calc_shot_box(times_t *times, pt) {*/

/*}*/
