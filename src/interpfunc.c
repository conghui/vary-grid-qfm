#include <math.h>
#include <stdlib.h>

#define MAX(a,b) ((a) < (b) ? (b) : (a))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/************************************************************************
*                          subroutine toep                              *
*************************************************************************
*                 toeplitz system solver: solves rf=g                   *
*************************************************************************
*/
static void toep (int m, const float *r,float *f,float *g,float *a){
  int i,j,jh;
  double c,e,v,w,bot;

  a[0]=1.;
  v=r[0];
  f[0]=g[0]/r[0];

  for (j=1; j<m; j++) {

    /* solve ra=v as in claerbout, fgdp, p. 57 */
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
*                          subroutine mksinc                            *
*************************************************************************
* derives tapered sinc interpolator coefficients by least squares       *
* spectral matching.  theory in wgc technical document by larner, 1979. *
*  dave hale, 1/31/83                                                   *
*************************************************************************
*   sinc    tapered sinc interpolation coefficients                     *
*   lsinc   length of sinc -- must be even; lsinc<20 recommended!       *
*   d       fractional distance to interpolated point; 0<=d<=1 required!*
*   space   workspace of lsinc*3 floats                                 *
*************************************************************************
*after calling mksinc, use the coefficients to interpolate as follows:  *
*                    j=lsinc-1                                          *
* y(t) = y(i+d) =      sum     sinc[j]*y[i+j-lsinc/2+1]                 *
*                      j=0                                              *
************************************************************************/
static void mksinc(float *sinc,int lsinc,float d,float *space)
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

  /* compute coefficients of toeplitz linear system */
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

void make_sinc_table(float *sinc_table, int ntab, int lsinc) {
  float dsinc=1./(ntab-1);
  float *space = (float *) malloc ( lsinc * 3 * sizeof(float) );

  /* contruct sinc table */
  for(int i=0; i< ntab; i++){
    float d=i*dsinc;
    mksinc( &sinc_table[i*lsinc], lsinc, d, space);
  }

  free(space);
}

static void sinc_interp_init(
    int *iaxis1, int *iaxis2, int *iaxis3, 
    int *spt1, int *spt2, int *spt3, 
    int ntab,
    float o1, float d1,
    float o2, float d2,
    float o3, float d3,
    int n1out, float o1out, float d1out,
    int n2out, float o2out, float d2out,
    int n3out, float o3out, float d3out)
{
  /*now lets allocate the tables*/
  float dsinc=1./(ntab-1);
  int i;
  float loc;

  /* do pointers for axis 1 */
  for(i=0; i< n1out; i++){
    loc=((o1out+d1out*i)-o1)/d1;
    iaxis1[i]=loc;
    spt1[i]=((loc-iaxis1[i])/dsinc)+.5;
  }

  /* do pointers for axis 2 */
  for(i=0; i< n2out; i++){
    loc=((o2out+d2out*i)-o2)/d2;
    iaxis2[i]=loc;
    spt2[i]=((loc-iaxis2[i])/(dsinc))+.5;
  }

  /* do pointers for axis 3 */
  for(i=0; i< n3out; i++){
    loc=((o3out+d3out*i)-o3)/d3;
    iaxis3[i]=loc;
    spt3[i]=((loc-iaxis3[i])/(dsinc))+.5;
  }

  /*first deal with possible boundary problem*/
  for(i=0;i<n1out;i++)if(spt1[i]>=ntab ||  spt1[i]<0) spt1[i]=0;
  for(i=0;i<n2out;i++)if(spt2[i]>=ntab ||  spt2[i]<0) spt2[i]=0;
  for(i=0;i<n3out;i++)if(spt3[i]>=ntab ||  spt3[i]<0) spt3[i]=0;
}

static void sinc_interp3d_(const float *input, float *output, const float *sinc_table,
    const int *iaxis1, const int *iaxis2, int *iaxis3, 
    const int *spt1, const int *spt2, const int *spt3, 
    int lsinc, int n1, int n2, int n3, int n1out, int n2out, int n3out)
{
#pragma omp parallel for schedule(dynamic, 1)
  for(int i3=0; i3 < n3out; i3++){
    for(int i2=0; i2 < n2out; i2++){
      for(int i1=0; i1 < n1out; i1++){
        float val=0.;
        int loc1,loc2,loc3;
        float f1,f2,f3;
        for(int c=0; c<lsinc; c++){
            loc3=iaxis3[i3]+c-lsinc/2+1;
            loc3=MIN(MAX(loc3,0),n3-1);
            f3=sinc_table[spt3[i3]*lsinc+c];
          for(int b=0; b<lsinc; b++){
              f2=sinc_table[spt2[i2]*lsinc+b];
              loc2=iaxis2[i2]+b-lsinc/2+1;
              loc2=MIN(MAX(loc2,0),n2-1);

            for(int a=0; a<lsinc; a++){
                f1=sinc_table[spt1[i1]*lsinc+a];
                loc1=iaxis1[i1]+a-lsinc/2+1;
                loc1=MIN(MAX(loc1,0),n1-1);
              val+=f1*f2*f3*input[loc1+loc2*n1+loc3*n2*n1];
            }
          }
        }
        output[i1+i2*n1out+i3*n1out*n2out]=val;
      }
    }
  }
}

void sinc_interp3d_1(const float *input, float *output, const float *sinc_table,
    int ntab, int lsinc,
    int n1, float o1, float d1,
    int n2, float o2, float d2,
    int n3, float o3, float d3,
    int n1out, float o1out, float d1out,
    int n2out, float o2out, float d2out,
    int n3out, float o3out, float d3out)
{

  int *iaxis1=(int*) malloc(n1out*sizeof(int));
  int *iaxis2=(int*) malloc(n2out*sizeof(int));
  int *iaxis3=(int*) malloc(n3out*sizeof(int));
  int *spt1=(int*) malloc(n1out*sizeof(int));
  int *spt2=(int*) malloc(n2out*sizeof(int));
  int *spt3=(int*) malloc(n3out*sizeof(int));

  sinc_interp_init(iaxis1,  iaxis2,  iaxis3, spt1,  spt2,  spt3, ntab, o1, d1, o2, d2, o3, d3, n1out, o1out, d1out, n2out, o2out, d2out, n3out, o3out, d3out);

  sinc_interp3d_( input, output,  sinc_table, iaxis1,  iaxis2, iaxis3, spt1,  spt2,  spt3, lsinc, n1, n2, n3, n1out, n2out, n3out);

  free(iaxis1); free(iaxis2); free(iaxis3);
  free(spt1); free(spt2); free(spt3);
}
