#include "resample.h"
#include "check.h"
#include <assert.h>

static int gs_ns;
static int gs_npts;
static float **gs_sinc_table;

/************************************************************************
*                          Subroutine toep                              *
*************************************************************************
*                 Toeplitz system solver: solves rf=g                   *
*************************************************************************
*/
static void toep (int m, const float *r,float *f,float *g,float *a){
  int i,j,jh;
  double c,e,v,w,bot;

  a[0]=1.;
  v=r[0];
  f[0]=g[0]/r[0];

  for (j=1; j<m; j++) {

    /* solve ra=v as in Claerbout, FGDP, p. 57 */
    e = a[j] = f[j] = 0.;
    for (i=0; i<j; i++)
    e += a[i]*r[j-i];
    c = e/v;
    v -= e*c;
    jh = j/2;
    for (i=0; i<=jh; i++) {
      bot = a[j-i]-c*a[i];
      a[i] -= c*a[j-i];
      a[j-i] = bot;
    }

    /* use a and v above to get f[i], i = 0,1,2,...,j */
    w = 0;
    for (i=0; i<j; i++)
      w += f[i]*r[j-i];
    c = (g[j]-w)/v;
    for (i=0; i<=j; i++)
      f[i] += c*a[j-i];
  }
}

/************************************************************************
*                          Subroutine mksinc                            *
*************************************************************************
* Derives tapered sinc interpolator coefficients by least squares       *
* spectral matching.  Theory in WGC technical document by Larner, 1979. *
*  Dave Hale, 1/31/83                                                   *
*************************************************************************
*   sinc    tapered sinc interpolation coefficients                     *
*   lsinc   length of sinc -- MUST BE EVEN; lsinc<20 RECOMMENDED!       *
*   d       fractional distance to interpolated point; 0<=d<=1 REQUIRED!*
*   space   workspace of lsinc*3 floats                                 *
*************************************************************************
*After calling mksinc, use the coefficients to interpolate as follows:  *
*                    j=lsinc-1                                          *
* y(t) = y(i+d) =      sum     sinc[j]*y[i+j-lsinc/2+1]                 *
*                      j=0                                              *
************************************************************************/
static void mksincit (float *sinc,int lsinc,float d,float *space)
{
  int j;
  float *b,*c,*work,pi,pi2,snyq,snat,s0,ds,eta,s;

  /* compute constants */
  pi = 3.141592654;
  pi2 = pi*2;
  snyq = 0.5;
  snat = snyq*(0.066+0.265*log((double)lsinc));
  s0 = 0.0;
  ds = (snat-s0)/(lsinc*2-1);
  eta = lsinc/2-1.0+d;

  /* segment work space */
  b = space; c = b+lsinc; work = c+lsinc;

  /* compute coefficients of Toeplitz linear system */
  for (j=0; j<lsinc; j++){
    for (s=s0,b[j]=c[j]=0.0; s<=snat; s+=ds)
      {
      b[j] += cos(pi2*s*j);
      c[j] += cos(pi2*s*(eta-j));
      }
  }

  /* solve the system for sinc coefficients */
  toep (lsinc,b,sinc,c,work);
}

static void mksinc_table(int npts, int nsinc, float **table)
{
  static float space[10000] = {0};

  float dtab = 1.0 / npts;
  float t = 0;
  for (int i = 0; i < npts; i++) {
    mksincit(table[i], nsinc, t, space);
    t = i * dtab;
  }
}


static void x_b_e(int oldn, float oldo, float oldd, int newn, float newo,  float newd, float *x, int *b, int *e)
{
  /*sf_warning("oldn, oldo, oldd: %d, %f, %f", oldn, oldo, oldd);*/
  /*sf_warning("newn, newo, newd: %d, %f, %f", newn, newo, newd);*/

  for (int i = 0; i < newn; i++) {
    x[i] = (newo + newd * i - oldo) / oldd + 1;
    /*sf_warning("x[%d]: %f", i, x[i]);*/
    if (x[i] > 0.999 && i < *b) {
      *b = i;
    }
    if (x[i] < oldn + 0.001 && i > *e) {
      *e = i;
    }

    /*break;*/
  }
}

static void calc_t(const float *x, int *k, int *t, int n, int sinc_npts)
{
  for (int i = 0; i < n; i++) {
    k[i] = floorf(x[i]);
    t[i] = roundf((x[i] - k[i]) * sinc_npts) + 0; // the value of t[i] will be used as index for sinc_table. So, counting from 0
  }
}

static void up_clip(int *t, int n, int value)
{
  for (int i = 0; i < n; i++) {
    t[i] = SF_MIN(t[i], value);
  }
}

static void updown_cip(int **t, int n1, int n2, int down, int up) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i2 = 0; i2 < n2; i2++) {
    for (int i1 = 0; i1 < n1; i1++) {
      if (t[i2][i1] > up) {
        t[i2][i1] = up;
      } else if (t[i2][i1] < down) {
        t[i2][i1] = down;
      }
    }
  }
}

void init_sinc_table(int nsinc, int npts)
{
  gs_ns         = nsinc;
  gs_npts       = npts;
  gs_sinc_table = sf_floatalloc2(nsinc, npts);

  mksinc_table(npts, nsinc, gs_sinc_table);

  /*sf_file ftable = sf_output("sinc.rsf");*/
  /*sf_putint(ftable, "n1", nsinc);*/
  /*sf_putint(ftable, "n2", npts);*/
  /*sf_putfloat(ftable, "o1", 0);*/
  /*sf_putfloat(ftable, "o2", 0);*/
  /*sf_putfloat(ftable, "d1", 0);*/
  /*sf_putfloat(ftable, "d2", 0);*/
  /*sf_floatwrite(gs_sinc_table[0], nsinc * npts, ftable);*/
}


void interpfield_(float ***oldf, float ***newf, bool extend,
    int on1, float oo1, float od1,  /* old */
    int on2, float oo2, float od2,
    int on3, float oo3, float od3,
    int nn1, float no1, float nd1,  /* new */
    int nn2, float no2, float nd2,
    int nn3, float no3, float nd3)
{
  float *x1 = sf_floatalloc(nn1);
  float *x2 = sf_floatalloc(nn2);
  float *x3 = sf_floatalloc(nn3);
  int b1 = nn1;
  int b2 = nn2;
  int b3 = nn3;
  int e1 = 0;
  int e2 = 0;
  int e3 = 0;
  int ns = gs_ns;

  memset(newf[0][0], 0, sizeof(float)*nn1*nn2*nn3);

  x_b_e(on1, oo1, od1, nn1, no1, nd1, x1, &b1, &e1);
  x_b_e(on2, oo2, od2, nn2, no2, nd2, x2, &b2, &e2);
  x_b_e(on3, oo3, od3, nn3, no3, nd3, x3, &b3, &e3);

  /*sf_warning("b1,b2,b3: %d, %d, %d", b1, b2, b3);*/
  /*sf_warning("e1,e2,e3: %d, %d, %d", e1, e2, e3);*/
  /*for (int i = 0; i < nn3; i++) {*/
    /*sf_warning("%f", x3[i]);*/
  /*}*/

  int *k1  = sf_intalloc(nn1); assert(k1);
  int *k2  = sf_intalloc(nn2);
  int *k3  = sf_intalloc(nn3);
  int *t1    = sf_intalloc(nn1);
  int *t2    = sf_intalloc(nn2);
  int *t3    = sf_intalloc(nn3);
  int **c1 = sf_intalloc2(ns, nn1);
  int **c2 = sf_intalloc2(ns, nn2);
  int **c3 = sf_intalloc2(ns, nn3);

  calc_t(x1, k1, t1, nn1, gs_npts);
  calc_t(x2, k2, t2, nn2, gs_npts);
  calc_t(x3, k3, t3, nn3, gs_npts);


  up_clip(t1, nn1, gs_npts - 1);
  up_clip(t2, nn2, gs_npts - 1);
  up_clip(t3, nn3, gs_npts - 1);

  /*sf_warning("sumk1, sumk2, sumk3: %d, %d, %d", isum1(k1,nn1), isum1(k2,nn2), isum1(k3, nn3));*/
  /*sf_warning("sumt1, sumt2, sumt3: %d, %d, %d", isum1(t1,nn1), isum1(t2,nn2), isum1(t3, nn3));*/

  for (int is = 0; is < ns; is++) {
    for (int i1 = 0; i1 < nn1; i1++) {
      c1[i1][is] = k1[i1] - 4 + is;
    }
    for (int i2 = 0; i2 < nn2; i2++) {
      c2[i2][is] = k2[i2] - 4 + is;
    }
    for (int i3 = 0; i3 < nn3; i3++) {
      c3[i3][is] = k3[i3] - 4 + is;
    }
  }

  updown_cip(c1, ns, nn1, 0, on1-1);
  updown_cip(c2, ns, nn2, 0, on2-1);
  updown_cip(c3, ns, nn3, 0, on3-1);

  /*for (int i = 0; i < ns; i++) {*/
    /*sf_warning("c1[0][%d]: %d", i, c1[0][i]);*/
  /*}*/
  /*sf_warning("write c1, c2, c3");*/
  /*write2di("c1.rsf", c1, ns, nn1);*/
  /*write2di("c2.rsf", c2, ns, nn2);*/
  /*write2di("c3.rsf", c3, ns, nn3);*/

  /*write1di("t1.rsf", t1, nn1);*/
  /*exit(0);*/

  if (extend) {
    /*sf_warning("computation of interpolation with extend = true");*/
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i3 = 0; i3 < nn3; i3++) {
      for (int i2 = 0; i2 < nn2; i2++) {
        for (int i1 = 0; i1 < nn1; i1++) {
          for (int ic = 0; ic < ns; ic++) {
            for (int ib = 0; ib < ns; ib++) {
              for (int ia = 0; ia < ns; ia++) {
                newf[i3][i2][i1] +=
                  oldf[c3[i3][ic]][c2[i2][ib]][c1[i1][ia]] *
                  gs_sinc_table[t1[i1]][ia] *
                  gs_sinc_table[t2[i2]][ib] *
                  gs_sinc_table[t3[i3]][ic];
              }
            }
          }
        }
      }
    }
  } else {
    /*sf_warning("computation of interpolation without extend");*/
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i3 = b3; i3 <= e3; i3++) {
      for (int i2 = b2; i2 <= e2; i2++) {
        for (int i1 = b1; i1 <= e1; i1++) {
          for (int ic = 0; ic < ns; ic++) {
            for (int ib = 0; ib < ns; ib++) {
              for (int ia = 0; ia < ns; ia++) {
                newf[i3][i2][i1] +=
                  oldf[c3[i3][ic]][c2[i2][ib]][c1[i1][ia]] *
                  gs_sinc_table[t1[i1]][ia] *
                  gs_sinc_table[t2[i2]][ib] *
                  gs_sinc_table[t3[i3]][ic];
              }
            }
          }
        }
      }
    }

  }
  free(x1); free(x2); free(x3);
  free(k1); free(k2); free(k3);
  free(t1); free(t2); free(t3);
  free(*c1); free(*c2); free(*c3);
  free(c1); free(c2); free(c3);
}

