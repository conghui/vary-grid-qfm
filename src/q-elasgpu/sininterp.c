#include <rsf.h>

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

int main(int argc, char **argv) {
  sf_init(argc, argv);

  /* information about input data set */
  sf_file Fin = sf_input("in");
  sf_file Fout = sf_output("out");

  int n1 = 1; int n2=1; int n3=1;
  float o1 = 0; float o2=0; float o3=0; 
  float d1 = 0; float d2=1; float d3=1; 

  int nn[10];
  float oo[10];
  float dd[10];
  for (int i = 1; i < 3; i++) {
    sf_axis ia = sf_iaxa(Fin, i); nn[i] = sf_n(ia); oo[i] = sf_o(ia); dd[i] = sf_d(ia); sf_raxa(ia);
  }
  n1 = nn[1]; o1 = oo[1]; d1 = dd[1];
  n2 = nn[2]; o2 = oo[2]; d2 = dd[2];
  n3 = nn[3]; o3 = oo[3]; d3 = dd[3];


  /* get interpolation parameters */
  int n1out, n2out, n3out;
  float o1out, o2out, o3out;
  float d1out, d2out, d3out;
  int lsinc, ntab;
  if(!sf_getfloat("o1out",&o1out)) o1out=o1;
  if(!sf_getfloat("d1out",&d1out)) d1out=d1;
  if(!sf_getfloat("o2out",&o2out)) o2out=o2;
  if(!sf_getfloat("d2out",&d2out)) d2out=d2;
  if(!sf_getfloat("o3out",&o3out)) o3out=o3;
  if(!sf_getfloat("d3out",&d3out)) d3out=d3;
  if(!sf_getint("n1out",&n1out)) n1out=((o1+d1*(n1-1)-o1out)/d1out+1);
  if(!sf_getint("n2out",&n2out)) n2out=((o2+d2*(n2-1)-o2out)/d2out+1);
  if(!sf_getint("n3out",&n3out)) n3out=((o3+d3*(n3-1)-o3out)/d3out+1);
  if(!sf_getint("ntab",&ntab)) ntab=101;
  if(!sf_getint("lsinc",&lsinc)) lsinc=10;

  sf_putint(Fout, "n1", n1out); sf_putfloat(Fout, "o1", o1out); sf_putfloat(Fout, "d1", d1out);
  sf_putint(Fout, "n2", n2out); sf_putfloat(Fout, "o2", o2out); sf_putfloat(Fout, "d2", d2out);
  sf_putint(Fout, "n3", n3out); sf_putfloat(Fout, "o3", o3out); sf_putfloat(Fout, "d3", d3out);

  int in_size = n1 * n2 * n3;
  int out_size = n1out * n2out * n3out;

  int *iaxis1=(int*) malloc(n1out*sizeof(int));
  int *iaxis2=(int*) malloc(n2out*sizeof(int));
  int *iaxis3=(int*) malloc(n3out*sizeof(int));

  int *spt1=(int*) malloc(n1out*sizeof(int));
  int *spt2=(int*) malloc(n2out*sizeof(int));
  int *spt3=(int*) malloc(n3out*sizeof(int));
  float *sinc_table=(float*) malloc(lsinc*sizeof(float)*ntab);
  float dsinc=1./(ntab-1);
  float nlen1=lsinc; 
  float nlen2=lsinc; 
  float nlen3=lsinc; 

  /* DO POINTERS FOR AXIS 1 */
  for(int i=0; i< n1out; i++){
    int loc=((o1out+d1out*i)-o1)/d1;
    iaxis1[i]=loc;
    spt1[i]=((loc-iaxis1[i])/dsinc)+.5;
  }

  /* DO POINTERS FOR AXIS 2 */
  for(int i=0; i< n2out; i++){
    int loc=((o2out+d2out*i)-o2)/d2;
    iaxis2[i]=loc;
    spt2[i]=((loc-iaxis2[i])/(dsinc))+.5;
  }

  /* DO POINTERS FOR AXIS 3 */
  for(int i=0; i< n3out; i++){
    int loc=((o3out+d3out*i)-o3)/d3;
    iaxis3[i]=loc;
    spt3[i]=((loc-iaxis3[i])/(dsinc))+.5;
  }

  float *space = (float *) malloc ( lsinc * 3 * sizeof(float) );
  /*first deal with possible boundary problem*/
  for(int i=0;i<n1out;i++)if(spt1[i]>=ntab ||  spt1[i]<0) spt1[i]=0;
  for(int i=0;i<n2out;i++)if(spt2[i]>=ntab ||  spt2[i]<0) spt2[i]=0;
  for(int i=0;i<n3out;i++)if(spt3[i]>=ntab ||  spt3[i]<0) spt3[i]=0;
  /* contruct sinc table */
  for(int i=0; i< ntab; i++){
    float d=i*dsinc;
    mksinc( &sinc_table[i*lsinc], lsinc, d, space);
  }

  float *input=(float*)malloc(in_size*sizeof(float));
  float *output=(float*)malloc(out_size*sizeof(float));

  for(int i3=0; i3 < n3out; i3++){
    for(int i2=0; i2 < n2out; i2++){
      for(int i1=0; i1 < n1out; i1++){
        float val=0.;
        int loc1,loc2,loc3;
        float f1,f2,f3;
        for(int c=0; c<nlen3; c++){
          if(nlen3==1){
            f3=1; loc3=i3;
          }
          else{
            loc3=iaxis3[i3]+c-lsinc/2+1;
            loc3=SF_MIN(SF_MAX(loc3,0),n3-1);
            f3=sinc_table[spt3[i3]*lsinc+c];
          }
          for(int b=0; b<nlen2; b++){
            if(nlen2==1){
              f2=1.;
              loc2=i2;
            } else{
              f2=sinc_table[spt2[i2]*lsinc+b];
              loc2=iaxis2[i2]+b-lsinc/2+1;
              loc2=SF_MIN(SF_MAX(loc2,0),n2-1);
            }

            for(int a=0; a<nlen1; a++){
              if(nlen1==1){
                f1=1.;
                loc1=i1;
              } else{
                f1=sinc_table[spt1[i1]*lsinc+a];
                loc1=iaxis1[i1]+a-lsinc/2+1;
                loc1=SF_MIN(SF_MAX(loc1,0),n1-1);
              }
              val+=f1*f2*f3*input[loc1+loc2*n1+loc3*n2*n1];
            }
          }
        }
        output[i1+i2*n1+i3*n1*n2]=val;
      }
    }
  }
  sf_floatwrite(output, out_size, Fout);

  return 0;
}
