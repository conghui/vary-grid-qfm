#include<math.h>
#include<stdio.h>
void mksinc (float *sinc,int lsinc,float d,float *space);;
void toep (int m,float *r,float *f,float *g,float *a);

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
void mksincit (float *sinc,int lsinc,float d,float *space)
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


/************************************************************************
*                          Subroutine toep                              *
*************************************************************************
*                 Toeplitz system solver: solves rf=g                   *
*************************************************************************
*/

void toep (int m,float *r,float *f,float *g,float *a){
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
