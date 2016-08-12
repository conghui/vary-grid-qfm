/*$

=head1 NAME

Interp - Interpolate dataset using sinc, linearm or nearest neighbor

=head1 SYNOPSIS

Interp < in.H o1 pars > out.H

=head1 INPUT PARAMETERS

=over 4

=item lsinc - integer

  [10]: length of interpolation operator
      (recommend < 20; MUST BE EVEN)

=item type  - int

   [2]: Type of interpolation,0 nearest neighbor,1-linear,2-sinc

=item o1out - float

  [o1]: First sample on axis1

=item o2out - float

  [o2]: First sample on axis2

=item o3out - float

  [o3]: First sample on axis3

=item d1out - float

  [d1]: Sampling of the output axis 1

=item d2out - float

  [d2]: Sampling of the output axis 2

=item d3out - float

  [d3]: Sampling of the output axis 3

=item n1out - int

  [max in/dout]: Number of samples in axis 1

=item n2out - int

  [max in/dout]: Number of samples in axis 2

=item n3out - int

  [max in/dout]: Number of samples in axis 3

=item maxsize-int

  [20]:  Amount of memory to use in megabytes

=item ntab  - int

  [101]: Interpolation table size (aka if outspace
      corresponds to inspace .012 .01 table will be chosen)

=back

=head1 DESCRIPTION

Interpolate dataset using sinc, linear, or nearest neighbor,
up to 3 dimensions(if it can be held in memory). If any of the
n1out,n2out,n3out, o1out,o2out,o3out or d1out,d2out,d3out is
omitted the corresponding value in the input data is used.

=head1 EXAMPLE

  Interp < in.H lsinc=12 type=2 > out.H
  conputes a 12-point sinc-interpolator on the input data.
  All of the standard n's, d's and o's are taken from the
  input data

=head1 CATEGORY

B<seis/filter>

=cut

>*/
/*
D. Rothman, 27Nov83
 * Keyword: sinc-interpolation
 *
 * Modified: May 26, 1998 - Bob - Changed to direct reference
 *                          to sinc rather than through psinc because
 *                          of 64Bit SGI
 * Modified: Aug 1, 1999   Bob Conbined with Interp2, extended options,
 *                          etc
 *
 */

#include <rsf.h>

#define YES 1
#define NO 0
#define NEAREST 0
#define LINEAR 1
#define SINC 2

static void putlin(const char *str) {
  sf_warning(str);
}

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


int main(int argc, char **argv)
{
  int lsinc, i;
  float *space;
  char temp_ch[256];
  int n1,n2,n3,n1out,n2out,n3out,nbig;
  float d,d1out,d2out,d3out,o1out,o2out,o3out,o1,o2,o3,d1,d2,d3;
  int type,do_axis[3],ntab,verb,maxsize;
  int *iaxis1,*iaxis2,*iaxis3;
  int *spt1 = NULL,*spt2=NULL,*spt3=NULL;
  float *faxis1=NULL,*faxis2=NULL,*faxis3=NULL;
  int nlen1,nlen2,nlen3;
  float *sinc_table=NULL,dsinc=0;
  float *input,*output;
  int in_size,out_size;
  int n1_lin,n2_lin,n3_lin;
  int n1_lout,n2_lout,n3_lout;
  int esize,nmax,unit;
  int iloc1,iloc2,iloc3;
  float loc;
  int i1,i2,i3;
  int b3,b2,b1,e3,e2,e1;

  sf_init(argc, argv);

  /* information about input data set */
  sf_file Fin = sf_input("in");
  sf_file Fout = sf_output("out");

  esize = sf_esize(Fin);
  if(esize!=4) sf_error("Can only deal with float data \n");

  nbig=1;
  o2=0;o3=0;d2=1;d3=1;n2=1;n3=1;

  int nn[10];
  float oo[10];
  float dd[10];
  for (int i = 1; i < 10; i++) {
    sf_axis ia = sf_iaxa(Fin, i); nn[i] = sf_n(ia); oo[i] = sf_o(ia); dd[i] = sf_d(ia); sf_raxa(ia);
    if (i > 3) {
      nbig *= nn[i];
    }
  }
  sf_warning("nbig: %d", nbig);
  n1 = nn[1]; o1 = oo[1]; d1 = dd[1];
  n2 = nn[2]; o2 = oo[2]; d2 = dd[2];
  n3 = nn[3]; o3 = oo[3]; d3 = dd[3];

  /*sf_axis ia1 = sf_iaxa(Fin, 1); n1 = sf_n(ia1); o1 = sf_o(ia1); d1 = sf_d(ia1); sf_raxa(ia1);*/
  /*sf_axis ia2 = sf_iaxa(Fin, 2); n2 = sf_n(ia2); o2 = sf_o(ia2); d2 = sf_d(ia2); sf_raxa(ia2);*/
  /*sf_axis ia3 = sf_iaxa(Fin, 3); n3 = sf_n(ia3); o3 = sf_o(ia3); d3 = sf_d(ia3); sf_raxa(ia3);*/
  /*sf_axis ia4 = sf_iaxa(Fin, 4);int n4 = sf_n(ia4);float o4 = sf_o(ia4);float d4 = sf_d(ia4); sf_raxa(ia4);*/
  /*sf_axis ia5 = sf_iaxa(Fin, 5);int n5 = sf_n(ia5);float o5 = sf_o(ia5);float d5 = sf_d(ia5); sf_raxa(ia5);*/
  /*sf_axis ia6 = sf_iaxa(Fin, 6);int n6 = sf_n(ia6);float o6 = sf_o(ia6);float d6 = sf_d(ia6); sf_raxa(ia6);*/
  /*sf_axis ia7 = sf_iaxa(Fin, 7);int n7 = sf_n(ia7);float o7 = sf_o(ia7);float d7 = sf_d(ia7); sf_raxa(ia7);*/
  /*sf_axis ia8 = sf_iaxa(Fin, 8);int n8 = sf_n(ia8);float o8 = sf_o(ia8);float d8 = sf_d(ia8); sf_raxa(ia8);*/
  /*sf_axis ia9 = sf_iaxa(Fin, 9);int n9 = sf_n(ia9);float o9 = sf_o(ia9);float d9 = sf_d(ia9); sf_raxa(ia9);*/

  /*nbig = n1 * n2 * n3 * n4 **/

  /* get interpolation parameters */
  if(0==sf_getint("type",&type)) type=SINC;
  if(type>2 || type<0) sf_error("Invalid type (0-nint 1-linear 2-sync)\n");
  do_axis[0]=NO; do_axis[1]=NO; do_axis[2]=NO;
  if(0==sf_getfloat("o1out",&o1out)) o1out=o1;else do_axis[0]=YES;
  if(0==sf_getfloat("d1out",&d1out)) d1out=d1;else do_axis[0]=YES;
  if(0==sf_getfloat("o2out",&o2out)) o2out=o2;else do_axis[1]=YES;
  if(0==sf_getfloat("d2out",&d2out)) d2out=d2;else do_axis[1]=YES;
  if(0==sf_getfloat("o3out",&o3out)) o3out=o3;else do_axis[2]=YES;
  if(0==sf_getfloat("d3out",&d3out)) d3out=d3;else do_axis[2]=YES;
  if(0==sf_getint("n1out",&n1out)) n1out=((o1+d1*(n1-1)-o1out)/d1out+1); else do_axis[0]=YES;
  if(0==sf_getint("n2out",&n2out)) n2out=((o2+d2*(n2-1)-o2out)/d2out+1); else do_axis[1]=YES;
  if(0==sf_getint("n3out",&n3out)) n3out=((o3+d3*(n3-1)-o3out)/d3out+1); else do_axis[2]=YES;
  if(0==sf_getint("lsinc",&lsinc)) lsinc=10;
  if(0==sf_getint("ntab",&ntab)) ntab=101;
  if(0==sf_getint("verb",&verb)) verb=0;
  if(0==sf_getint("maxsize",&maxsize)) maxsize=2000;
  maxsize=maxsize*1000000;

  sf_putint(Fout, "n1", n1out); sf_putfloat(Fout, "o1", o1out); sf_putfloat(Fout, "d1", d1out);
  sf_putint(Fout, "n2", n2out); sf_putfloat(Fout, "o2", o2out); sf_putfloat(Fout, "d2", d2out);
  sf_putint(Fout, "n3", n3out); sf_putfloat(Fout, "o3", o3out); sf_putfloat(Fout, "d3", d3out);

  /*take care of the output space */
  if(do_axis[0]==NO && do_axis[1]==NO && do_axis[2]==NO)
    sf_error("You aren't interpolating along any axis \n");
  if(type==NEAREST) strcpy(temp_ch,"nearest neighbor");
  else if(type==LINEAR) strcpy(temp_ch,"linear interpolation");
  else sprintf(temp_ch,"sinc interpolation with %d points",lsinc);
  putlin(temp_ch);
  if(verb==1) fprintf(stderr,"%s\n",temp_ch);
  if(do_axis[0]==YES){
    sprintf(temp_ch,"n1=%d o1=%f d1=%f",n1out,o1out,d1out);
    putlin(temp_ch);if(verb==1) fprintf(stderr,"%s\n",temp_ch);
  }
  if(do_axis[1]==YES){
    sprintf(temp_ch,"n2=%d o2=%f d2=%f",n2out,o2out,d2out);
    putlin(temp_ch);if(verb==1) fprintf(stderr,"%s\n",temp_ch);
  }
  if(do_axis[2]==YES){
    sprintf(temp_ch,"n3=%d o3=%f d3=%f",n3out,o3out,d3out);
    putlin(temp_ch);if(verb==1) fprintf(stderr,"%s\n",temp_ch);
  }

  if(do_axis[2]==NO){
    nbig=n3*nbig;
    nmax=2;
    if(do_axis[1]==NO){
      nbig=nbig*n2; nmax=1; in_size=n1;out_size=n1out;
      n2_lin=1; n3_lin=1;
      n2_lout=1; n3_lout=1;
    }
    else{
      in_size=(n1*n2);out_size=(n1out*n2out);
      n3_lin=1; n2_lin=n2;
      n3_lout=1; n2_lout=n2out;
    }
  }
  else{
    in_size=(n1*n2*n3);out_size=n1out*n2out*n3out;
    n3_lout=n3out; n2_lout=n2out;
    n3_lin=n3; n2_lin=n2;
  }

  n1_lout=n1out;
  n1_lin=n1;

  unit=in_size+out_size;
  if(unit*esize  > maxsize)
    sf_error("maxsize is not big enough, must be %d with current pars\n",
        unit*esize/1000/1000);

  /*NOW LETS ALLOCATE THE TABLES*/
  iaxis1=(int*) malloc(n1out*sizeof(int));
  iaxis2=(int*) malloc(n2out*sizeof(int));
  iaxis3=(int*) malloc(n3out*sizeof(int));
  if(type==LINEAR){
    faxis1=(float*) malloc(n1out*sizeof(float));
    faxis2=(float*) malloc(n2out*sizeof(float));
    faxis3=(float*) malloc(n3out*sizeof(float));
    if(do_axis[0]==YES)nlen1=2; else nlen1=1;
    if(do_axis[1]==YES)nlen2=2; else nlen2=1;
    if(do_axis[2]==YES)nlen3=2; else nlen3=1;
  }
  else if(type==SINC){
    spt1=(int*) malloc(n1out*sizeof(int));
    spt2=(int*) malloc(n2out*sizeof(int));
    spt3=(int*) malloc(n3out*sizeof(int));
    sinc_table=(float*) malloc(lsinc*sizeof(float)*ntab);
    dsinc=1./(ntab-1);
    if(do_axis[0]==YES)nlen1=lsinc; else nlen1=1;
    if(do_axis[1]==YES)nlen2=lsinc; else nlen2=1;
    if(do_axis[2]==YES)nlen3=lsinc; else nlen3=1;
  }
  else{ nlen1=1; nlen2=1; nlen3=1;}


  /* DO POINTERS FOR AXIS 1 */
  for(i=0; i< n1out; i++){
    loc=((o1out+d1out*i)-o1)/d1;
    iaxis1[i]=loc;
    if(do_axis[0]==YES){
      if(type==LINEAR) faxis1[i]=loc-iaxis1[i];
      else if(type==SINC) spt1[i]=((loc-iaxis1[i])/dsinc)+.5;
    }
    else{
      if(type==LINEAR) faxis1[i]=0.;
      else if(type==SINC) spt1[i]=1;
    }
  }

  /* DO POINTERS FOR AXIS 2 */
  for(i=0; i< n2out; i++){
    loc=((o2out+d2out*i)-o2)/d2;
    iaxis2[i]=loc;
    if(do_axis[1]==YES){
      if(type==LINEAR) faxis2[i]=loc-iaxis2[i];
      else if(type==SINC) spt2[i]=((loc-iaxis2[i])/(dsinc))+.5;
    }
    else{
      if(type==LINEAR) faxis2[i]=0.;
      else if(type==SINC) spt2[i]=1;
    }
  }

  /* DO POINTERS FOR AXIS 3 */
  for(i=0; i< n3out; i++){
    loc=((o3out+d3out*i)-o3)/d3;
    iaxis3[i]=loc;
    if(do_axis[2]==YES){
      if(type==LINEAR) faxis3[i]=loc-iaxis3[i];
      else if(type==SINC) spt3[i]=((loc-iaxis3[i])/(dsinc))+.5;
    }
    else{
      if(type==LINEAR) faxis3[i]=0.;
      else if(type==SINC) spt3[i]=1;
    }
  }

  if(verb) fprintf(stderr,"Finished constructing pointers \n");

  if(type==SINC){
    space = (float *) malloc ( lsinc * 3 * sizeof(float) );
    /*first deal with possible boundary problem*/
    for(i=0;i<n1out;i++)if(spt1[i]>=ntab ||  spt1[i]<0) spt1[i]=0;
    for(i=0;i<n2out;i++)if(spt2[i]>=ntab ||  spt2[i]<0) spt2[i]=0;
    for(i=0;i<n3out;i++)if(spt3[i]>=ntab ||  spt3[i]<0) spt3[i]=0;
    /* contruct sinc table */
    for(i=0; i< ntab; i++){
      d=i*dsinc;
      mksinc( &sinc_table[i*lsinc], lsinc, d, space);
    }
    if(verb)
      fprintf(stderr,"finished constructing sinc table of size %d by %d\n",
          lsinc,ntab);
  }

  /*now it is time to get to work doing interpolation */
  input=(float*)malloc(in_size*sizeof(float));
  output=(float*)malloc(out_size*sizeof(float));

  for (int ih = 0; ih < nbig; ih++) {
    sf_floatread(input, in_size, Fin);

    /*NEAREST NEIGHBOR CASE */
    if(type==NEAREST){
      for(i3=0; i3 < n3_lout; i3++){
        iloc3=SF_MIN(SF_MAX(iaxis3[i3]+1.5,0),n3_lin-1);
        for(i2=0; i2 < n2_lout; i2++){
          iloc2=SF_MIN(SF_MAX(iaxis2[i2]+1.5,0),n2_lin-1);
          for(i1=0; i1 < n1_lout; i1++){
            iloc1=SF_MIN(SF_MAX(iaxis1[i1]+1.5,0),n1_lin-1);
            output[i1+n1_lout*i2+i3*n1_lout*n2_lout]=
              input[iloc1+n1_lin*iloc2+iloc3*n1_lin*n2_lin];
          }
        }
      }
    }
    else if(type==LINEAR){
      for(i3=0; i3 < n3_lout; i3++){
        b3=SF_MAX(SF_MIN(iaxis3[i3],n3_lin-1),0);
        e3=SF_MIN(SF_MAX(iaxis3[i3]+1,0),n3_lin-1);
        for(i2=0; i2 < n2_lout; i2++){
          b2=SF_MAX(SF_MIN(iaxis2[i2],n2_lin-1),0);
          e2=SF_MIN(SF_MAX(iaxis2[i2]+1,0),n2_lin-1);
          for(i1=0; i1 < n1_lout; i1++){
            b1=SF_MAX(SF_MIN(iaxis1[i1],n1_lin-1),0);
            e1=SF_MIN(SF_MAX(iaxis1[i1]+1,0),n1_lin-1);
            output[i1+n1_lout*i2+i3*n1_lout*n2_lout]=
              (1.-faxis1[i1])*(1-faxis2[i2])*(1.-faxis3[i3])*
              input[b1+b2*n1_lin+b3*n1_lin*n2_lin]+
              (faxis1[i1])*(1-faxis2[i2])*(1.-faxis3[i3])*
              input[e1+b2*n1_lin+b3*n1_lin*n2_lin]+
              (1.-faxis1[i1])*(faxis2[i2])*(1.-faxis3[i3])*
              input[b1+e2*n1_lin+b3*n1_lin*n2_lin]+
              (faxis1[i1])*(faxis2[i2])*(1.-faxis3[i3])*
              input[e1+e2*n1_lin+b3*n1_lin*n2_lin]+
              (1.-faxis1[i1])*(1.-faxis2[i2])*(faxis3[i3])*
              input[b1+b2*n1_lin+e3*n1_lin*n2_lin]+
              (faxis1[i1])*(1.-faxis2[i2])*(faxis3[i3])*
              input[e1+b2*n1_lin+e3*n1_lin*n2_lin]+
              (1.-faxis1[i1])*(faxis2[i2])*(faxis3[i3])*
              input[b1+e2*n1_lin+e3*n1_lin*n2_lin]+
              (faxis1[i1])*(faxis2[i2])*(faxis3[i3])*
              input[e1+e2*n1_lin+e3*n1_lin*n2_lin];
          }
        }
      }
    }
    else{  /* sinc interpolation */
#pragma omp parallel for schedule(dynamic, 1)
      for(int i3=0; i3 < n3_lout; i3++){
        for(int i2=0; i2 < n2_lout; i2++){
          for(int i1=0; i1 < n1_lout; i1++){
            float val=0.;
            int loc1,loc2,loc3;
            float f1,f2,f3;
            for(int c=0; c<nlen3; c++){
              if(nlen3==1){
                f3=1; loc3=i3;
              }
              else{
                loc3=iaxis3[i3]+c-lsinc/2+1;
                loc3=SF_MIN(SF_MAX(loc3,0),n3_lin-1);
                f3=sinc_table[spt3[i3]*lsinc+c];
              }
              for(int b=0; b<nlen2; b++){
                if(nlen2==1){
                  f2=1.;
                  loc2=i2;
                } else{
                  f2=sinc_table[spt2[i2]*lsinc+b];
                  loc2=iaxis2[i2]+b-lsinc/2+1;
                  loc2=SF_MIN(SF_MAX(loc2,0),n2_lin-1);
                }

                for(int a=0; a<nlen1; a++){
                  if(nlen1==1){
                    f1=1.;
                    loc1=i1;
                  } else{
                    f1=sinc_table[spt1[i1]*lsinc+a];
                    loc1=iaxis1[i1]+a-lsinc/2+1;
                    loc1=SF_MIN(SF_MAX(loc1,0),n1_lin-1);
                  }
                  val+=f1*f2*f3*input[loc1+loc2*n1_lin+loc3*n2_lin*n1_lin];
                }
              }
            }
            output[i1+i2*n1_lout+i3*n1_lout*n2_lout]=val;
          }
        }
      }
    }

    sf_floatwrite(output, out_size, Fout);
  }


  return(0);
}

