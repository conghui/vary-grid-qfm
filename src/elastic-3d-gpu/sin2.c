
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
#include "interpfunc.h"

int main(int argc, char **argv)
{
  int lsinc;
  int n1,n2,n3,n1out,n2out,n3out;
  float d1out,d2out,d3out,o1out,o2out,o3out,o1,o2,o3,d1,d2,d3;
  int ntab;
  int *iaxis1,*iaxis2,*iaxis3;
  int *spt1 = NULL,*spt2=NULL,*spt3=NULL;
  float *sinc_table=NULL;
  float *input,*output;
  int in_size,out_size;

  sf_init(argc, argv);

  /* information about input data set */
  sf_file Fin = sf_input("in");
  sf_file Fout = sf_output("out");

  o2=0;o3=0;d2=1;d3=1;n2=1;n3=1;

  int nn[10];
  float oo[10];
  float dd[10];
  int nbig = 1;
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

  if(0==sf_getfloat("o1out",&o1out)) o1out=o1;
  if(0==sf_getfloat("d1out",&d1out)) d1out=d1;
  if(0==sf_getfloat("o2out",&o2out)) o2out=o2;
  if(0==sf_getfloat("d2out",&d2out)) d2out=d2;
  if(0==sf_getfloat("o3out",&o3out)) o3out=o3;
  if(0==sf_getfloat("d3out",&d3out)) d3out=d3;
  if(0==sf_getint("n1out",&n1out)) n1out=((o1+d1*(n1-1)-o1out)/d1out+1); 
  if(0==sf_getint("n2out",&n2out)) n2out=((o2+d2*(n2-1)-o2out)/d2out+1); 
  if(0==sf_getint("n3out",&n3out)) n3out=((o3+d3*(n3-1)-o3out)/d3out+1); 
  if(0==sf_getint("lsinc",&lsinc)) lsinc=8;
  if(0==sf_getint("ntab",&ntab)) ntab=1000;

  sf_putint(Fout, "n1", n1out); sf_putfloat(Fout, "o1", o1out); sf_putfloat(Fout, "d1", d1out);
  sf_putint(Fout, "n2", n2out); sf_putfloat(Fout, "o2", o2out); sf_putfloat(Fout, "d2", d2out);
  sf_putint(Fout, "n3", n3out); sf_putfloat(Fout, "o3", o3out); sf_putfloat(Fout, "d3", d3out);

  in_size=(n1*n2*n3);out_size=n1out*n2out*n3out;

  iaxis1=(int*) malloc(n1out*sizeof(int));
  iaxis2=(int*) malloc(n2out*sizeof(int));
  iaxis3=(int*) malloc(n3out*sizeof(int));
  spt1=(int*) malloc(n1out*sizeof(int));
  spt2=(int*) malloc(n2out*sizeof(int));
  spt3=(int*) malloc(n3out*sizeof(int));
  sinc_table=(float*) malloc(lsinc*sizeof(float)*ntab);

  /*now it is time to get to work doing interpolation */
  input=(float*)malloc(in_size*sizeof(float));
  output=(float*)malloc(out_size*sizeof(float));
  make_sinc_table(sinc_table, ntab, lsinc);

  for (int ib = 0; ib < nbig; ib++) {
    sf_floatread(input, in_size, Fin);
    sinc_interp3d_1(input, output, sinc_table, ntab, lsinc, n1, o1, d1, n2, o2, d2, n3, o3, d3, n1out, o1out, d1out, n2out, o2out, d2out, n3out, o3out, d3out);
    sf_floatwrite(output, out_size, Fout);
  }

  return(0);
}

