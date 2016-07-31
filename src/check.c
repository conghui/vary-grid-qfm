#include "check.h"
#include <rsf.h>

float fsum1(const float *a, int n) {
  float s = 0;
  for (int i = 0; i < n; i++) {
    s += a[i];
  }

  return s;
}

int isum1(const int *a, int n) {
  int s = 0;
  for (int i = 0; i < n; i++) {
    s += a[i];
  }

  return s;

}

void write2df(const char *file, float **dat, int n1, int n2) {
  sf_file f = sf_output(file);
  sf_putint(f, "n1", n1); sf_putfloat(f, "o1", 0); sf_putfloat(f, "d1", 0);
  sf_putint(f, "n2", n2); sf_putfloat(f, "o2", 0); sf_putfloat(f, "d2", 0);
  sf_floatwrite(dat[0], n1*n2, f);
}

void write2di(const char *file, int **dat, int n1, int n2) {
  sf_file f = sf_output(file);
  sf_putint(f, "n1", n1); sf_putint(f, "o1", 0); sf_putint(f, "d1", 0);
  sf_putint(f, "n2", n2); sf_putint(f, "o2", 0); sf_putint(f, "d2", 0);
  sf_putint(f, "n3", 1); 
  sf_settype(f, SF_INT);
  sf_intwrite(dat[0], n1*n2, f);
}

void write1di(const char *file, int *dat, int n1) {
  sf_file f = sf_output(file);
  sf_putint(f, "n1", n1); sf_putint(f, "o1", 0); sf_putint(f, "d1", 0);
  sf_putint(f, "n2", 1); 
  sf_settype(f, SF_INT);
  sf_intwrite(dat, n1, f);
}
