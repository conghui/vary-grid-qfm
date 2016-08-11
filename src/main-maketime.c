
#include <rsf.h>

#include "fastmarch.h"

int main(int argc, char* argv[])
{
  sf_init(argc,argv);

  int j2, j3;
  bool plane[3] = {false, false, false};
  const int nflag = 50000;
  int *flag = sf_intalloc(nflag);
  memset(flag, 0, sizeof *flag * nflag);
  if (!sf_getint("j2", &j2)) j2 = 1;
  if (!sf_getint("j3", &j3)) j3 = 1;

  sf_file Fvel = sf_input("in");
  sf_file Ftime = sf_output("out");

  sf_axis az = sf_iaxa(Fvel, 1); sf_setlabel(az, "z"); sf_raxa(az);
  sf_axis ax = sf_iaxa(Fvel, 2); sf_setlabel(ax, "x"); sf_raxa(ax);
  sf_axis ay = sf_iaxa(Fvel, 3); sf_setlabel(ay, "y"); sf_raxa(ay);

  int n1 = sf_n(az); float o1 = sf_o(az); float d1 = sf_d(az);
  int n2 = sf_n(ax); float o2 = sf_o(ax); float d2 = sf_d(ax);
  int n3 = sf_n(ay); float o3 = sf_o(ay); float d3 = sf_d(ay);

  int n4out = ceilf((float)n2 / (float)(j2));
  int n5out = ceilf((float)n3 / (float)(j3));
  sf_axis asx = sf_maxa(n4out, 0, d2*j2);
  sf_axis asy = sf_maxa(n5out, 0, d3*j3);

  sf_oaxa(Ftime, az, 1);
  sf_oaxa(Ftime, ax, 2);
  sf_oaxa(Ftime, ay, 3);
  sf_oaxa(Ftime, asx, 4);
  sf_oaxa(Ftime, asy, 5);

  float ***v = sf_floatalloc3(n1, n2, n3);
  sf_floatread(v[0][0], n1*n2*n3, Fvel);

  for (int i3 = 0; i3 < n3; i3++) {
    for (int i2 = 0; i2 < n2; i2++) {
      for (int i1 = 0; i1 < n1; i1++) {
        v[i3][i2][i1] = (1.0 / v[i3][i2][i1]) * (1.0 / v[i3][i2][i1]);
      }
    }
  }

  float *time = sf_floatalloc(n1 * n2 * n3);
  fastmarch_init(n3, n2, n1);

  int order = 3;

  for (int i3 = 0; i3 < n3; i3 += j3) {
    for (int i2 = 0; i2 < n2; i2 += j2) {
      float s1 = 0;
      float s2 = i2 * d2;
      float s3 = i3 * d3;

      sf_warning("i2: (%d/%d), i3: (%d/%d)", i2 % j2, n4out, i3 % j3, n5out);
      fastmarch(time,v[0][0],flag,plane,
          n3,n2,n1,o3,o2,o1,d3,d2,d1,
          s3,s2,s1,1,1,1,order);

      sf_floatwrite(time, n1*n2*n3, Ftime);
    }

  }

  free(time);
  free(**v); free(*v); free(v);

  fastmarch_close();

  exit(0);
}
