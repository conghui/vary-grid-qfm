/* 3D acoustic time-domain FD modeling
   2Nth order in space, 2nd order in time
   with optimized FD scheme option and hybrid one-way ABC option
   adj flag indicates backwards source injection, not exact adjoint propagator
   */
/*
   Copyright (C) 2014 Colorado School of Mines
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
   */
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
#endif
#include "fdutil.h"
#include "step-forward.h"
#include "common.h"

#ifndef SF_HAS_SSE
#undef __SSE__
#endif
#ifndef SF_HAS_AVX
#undef __AVX__
#endif

#ifdef sun
#define restrict
#endif

#if defined __SSE__ || defined __AVX__
#include <immintrin.h>
#endif


  int
main(int argc, char** argv)
{

  bool verb, fsrf, snap, expl, dabc, cden, adj;
  bool optfd, hybrid, sinc;
  int jsnap, jdata;

  /* I/O files */
  sf_file file_wav=NULL; /* wavelet */
  sf_file file_vel=NULL; /* velocity */
  sf_file file_den=NULL; /* density */
  sf_file file_wfl=NULL; /* wavefield */
  sf_file file_dat=NULL; /* data */
  sf_file file_src=NULL; /* sources */
  sf_file file_rec=NULL; /* receivers */

  /* cube axes */
  sf_axis at = NULL, az = NULL, ax = NULL, ay = NULL;
  sf_axis as = NULL, ar = NULL;

  int nbd;  /* ABC boundary size */
  int fdorder;  /* finite difference spatial accuracy order */
  int nzpad,nxpad,nypad; /* boundary padded model size */
  int ix,iy,it,nx,ny,nz,nt,ns,nr;
  float dx,dy,dz,dt,dt2;
  float* damp=NULL; /* damping profile for hybrid bc */
  float* ws;  /* wavelet */
  float*** vel=NULL;  /* velocity */
  float*** rho=NULL; /* density */
  float*** u0=NULL;  /* wavefield array u@t-1 (u@t+1) */
  float*** u1=NULL;  /* wavefield array u@t */
  float* u_dat=NULL; /* output data */
  float*** ptr_tmp=NULL;
  pt3d* src3d=NULL;  /* source position */
  pt3d* rec3d=NULL;  /*receiver position*/
  scoef3d cssinc = NULL, crsinc = NULL;
  lint3d cslint = NULL, crlint = NULL;

  /* FDM structure */
  fdm3d fdm = NULL;
  abcone3d abc = NULL;
  sponge spo = NULL;

  float* fdcoef_d2;
  float* fdcoef_d1;

  sf_axis acz = NULL, acx = NULL, acy = NULL;
  int nqz, nqx, nqy;
  float oqz, oqx, oqy, dqz, dqx, dqy;

  float** oslice = NULL; /* output 3D wavefield slice-by-slice */
  float*** tmp_array;

  double wall_clock_time_s, wall_clock_time_e;

  const int SECOND_DERIV = 2;
  const int FIRST_DERIV = 1;

  int nop;

#if defined _OPENMP && _DEBUG
  double tic;
  double toc;
#endif

  /* init RSF */
  sf_init(argc,argv);

#ifdef _OPENMP
  omp_init();
  wall_clock_time_s = omp_get_wtime();
#else
  wall_clock_time_s = (double) clock() / CLOCKS_PER_SEC;
#endif

  if (!sf_getbool("verb",&verb))  verb=false; /* Verbosity flag */
  if (!sf_getbool("snap",&snap))  snap=false; /* Wavefield snapshots flag */
  if (!sf_getbool("expl",&expl))  expl=false; /* Multiple sources, one wvlt*/
  if (!sf_getbool("dabc",&dabc))  dabc=false; /* Absorbing BC */
  if (!sf_getbool("cden",&cden))  cden=false; /* Constant density */
  if (!sf_getbool("adj",&adj))    adj=false; /* adjoint flag */

  if (!sf_getbool("free",&fsrf) && !sf_getbool("fsrf",&fsrf)) fsrf=false; /* Free surface flag */
  if (!sf_getbool("optfd",&optfd))  optfd=false; /* optimized FD coefficients flag */
  if (!sf_getint("fdorder",&fdorder))  fdorder=4; /* spatial FD order */
  if (!sf_getbool("hybridbc",&hybrid))  hybrid=false;  /* hybrid Absorbing BC */
  if (!sf_getbool("sinc",&sinc)) sinc=false; /* sinc source injection */

  /* Initialize variables */
  file_wav = sf_input("in"); /* wavelet */
  file_vel = sf_input("vel"); /* velocity */
  file_src = sf_input("sou"); /* sources */
  file_rec = sf_input("rec"); /* receivers */
  file_dat = sf_output("out"); /* data */

  if (snap)  file_wfl = sf_output("wfl"); /* wavefield */
  if (!cden) {
    if (sf_getstring("den")) {
      file_den = sf_input ("den"); /* density */
    } else {
      cden = true;
      if (verb) sf_warning("No density file provided, running with constant density");
    }
  }

  at = sf_iaxa(file_wav,2); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
  az = sf_iaxa(file_vel,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
  ax = sf_iaxa(file_vel,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */
  ay = sf_iaxa(file_vel,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /* space */

    if(verb && fsrf) sf_warning("free surface condition");

  sf_axis asz, asx, asy;
  asz = sf_iaxa(file_src, 2); asx = sf_iaxa(file_src, 3); asy = sf_iaxa(file_src, 4);
  as = sf_maxa(sf_n(asz) * sf_n(asx) * sf_n(asy), sf_d(asx), sf_o(asx));
  /*as = sf_iaxa(file_src,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); [> sources <]*/

  sf_axis arz, arx, ary;
  arz = sf_iaxa(file_rec, 2); arx = sf_iaxa(file_rec, 3); ary = sf_iaxa(file_rec, 4);
  ar = sf_maxa(sf_n(arz) * sf_n(arx) * sf_n(ary), sf_d(arx), sf_o(arx));
  /*ar = sf_iaxa(file_rec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); [> receivers <]*/

  nt = sf_n(at); dt = sf_d(at);
  nz = sf_n(az); dz = sf_d(az);
  nx = sf_n(ax); dx = sf_d(ax);
  ny = sf_n(ay); dy = sf_d(ay);
  ns = sf_n(as);
  nr = sf_n(ar);

  /* other execution parameters */
  if (snap) {
    if (!sf_getint("jsnap",&jsnap))  jsnap=nt;
    /* # of t steps at which to save wavefield */
  }
  if (!sf_getint("jdata",&jdata)) jdata=1;
  /* # of t steps at which to save receiver data */

  /* setup output data header */
  sf_oaxa(file_dat,ar,1);
  sf_setn(at,(nt-1)/jdata+1);
  sf_setd(at,dt*jdata);
  sf_oaxa(file_dat,at,2);

  /* wavefield cut params */
  /* setup output wavefield header */
  if (snap) {
    if (!sf_getint  ("nqz",&nqz)) nqz=sf_n(az); /* Saved wfld window nz */
    if (!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax); /* Saved wfld window nx */
    if (!sf_getint  ("nqy",&nqy)) nqy=sf_n(ay); /* Saved wfld window ny */

    if (!sf_getfloat("oqz",&oqz)) oqz=sf_o(az); /* Saved wfld window oz */
    if (!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax); /* Saved wfld window ox */
    if (!sf_getfloat("oqy",&oqy)) oqy=sf_o(ay); /* Saved wfld window oy */

    if (!sf_getfloat("dqz",&dqz)) dqz=sf_d(az); /* Saved wfld window dz */
    if (!sf_getfloat("dqx",&dqx)) dqx=sf_d(ax); /* Saved wfld window dx */
    if (!sf_getfloat("dqy",&dqy)) dqy=sf_d(ay); /* Saved wfld window dy */

    acz = sf_maxa(nqz,oqz,dqz); if (verb) sf_raxa(acz);
    acx = sf_maxa(nqx,oqx,dqx); if (verb) sf_raxa(acx);
    acy = sf_maxa(nqy,oqy,dqy); if (verb) sf_raxa(acy);
    /* check if the imaging window fits in the wavefield domain */
    sf_setn(at,(nt-1)/jsnap+1);
    sf_setd(at,dt*jsnap);
    if (verb) sf_raxa(at);

    sf_oaxa(file_wfl,acz,1);
    sf_oaxa(file_wfl,acx,2);
    sf_oaxa(file_wfl,acy,3);
    sf_oaxa(file_wfl,at,4);
  }

  /* 2-2N finite difference coefficient */
  nop = fdorder/2; /* fd half-length stencil */
  if (!sf_getint("nb",&nbd) || nbd<nop)  nbd=nop;
  if (dabc && hybrid && nbd<=nop) nbd = 2*nop;

  /* expand domain for FD operators and ABC */
  fdm = fdutil3d_init(verb,fsrf,az,ax,ay,nbd,1);

  sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if (verb) sf_raxa(az);
  sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if (verb) sf_raxa(ax);
  sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); if (verb) sf_raxa(ay);

  /* Precompute coefficients */
  dt2 = dt*dt;
  nzpad = nz+2*nbd;  nxpad = nx+2*nbd;  nypad = ny+2*nbd;

  fdcoef_d2 = compute_fdcoef(nop,dz,dx,dy,optfd,SECOND_DERIV);
  fdcoef_d1 = compute_fdcoef(nop,dz,dx,dy,optfd,FIRST_DERIV);

  /* Allocate memories */
  if (expl) ws = sf_floatalloc(1);
  else      ws = sf_floatalloc(ns);
  vel = sf_floatalloc3(nzpad,nxpad,nypad);
  if (!cden) rho = sf_floatalloc3(nzpad,nxpad,nypad);
  u_dat = sf_floatalloc(nr);
  src3d = pt3dalloc1(ns);
  rec3d = pt3dalloc1(nr);
  if (snap) oslice = sf_floatalloc2(sf_n(acz),sf_n(acx));

  /* source and receiver position */
  pt3dread1(file_src,src3d,ns,3);  /* read format: (x,y,z) */
  if (sinc) cssinc = sinc3d_make(ns,src3d,fdm);
  else      cslint = lint3d_make(ns,src3d,fdm);

  pt3dread1(file_rec,rec3d,nr,3);  /*read format: (x,y,z) */
  if (sinc) crsinc = sinc3d_make(nr,rec3d,fdm);
  else      crlint = lint3d_make(nr,rec3d,fdm);

  /* temperary array */
  tmp_array = sf_floatalloc3(nz,nx,ny);

  /* read velocity and pad */
  sf_floatread(tmp_array[0][0],nz*nx*ny,file_vel);
  expand3d(tmp_array,vel,fdm);
  /* read density and pad */
  if (!cden) {
    sf_floatread(tmp_array[0][0],nz*nx*ny,file_den);
    expand3d(tmp_array,rho,fdm);
  }

  free(**tmp_array);  free(*tmp_array);  free(tmp_array);

  /* A1 one-way ABC implicit scheme coefficients  */
  if (dabc) {
    abc = abcone3d_make(nbd,dt,vel,fsrf,fdm);
    if (hybrid)
      damp = damp_make(nbd-nop); /* compute damping profiles for hybrid bc */
    else
      spo = sponge_make(fdm->nb);
  }

  /* allocate memory for wavefield variables */
  u0 = sf_floatalloc3(nzpad,nxpad,nypad);
  u1 = sf_floatalloc3(nzpad,nxpad,nypad);

  /* initialize variables */
  memset(u0[0][0],0,sizeof(float)*nzpad*nxpad*nypad);
  memset(u1[0][0],0,sizeof(float)*nzpad*nxpad*nypad);
  memset(u_dat,0,sizeof(float)*nr);

  /* v = (v*dt)^2 */
  for (ix=0;ix<nzpad*nxpad*nypad;ix++)
    *(vel[0][0]+ix) *= *(vel[0][0]+ix)*dt2;
  if (fsrf && !hybrid) {
    for (iy=0; iy<nypad; iy++)
      for (ix=0; ix<nxpad; ix++)
        memset(vel[iy][ix],0,sizeof(float)*(fdm->nb+1));
  }

  sf_warning("ns: %d, nr: %d", ns, nr);
  sf_warning("run only one shot");
  for (it=0; it<nt; it++) {
    if (verb)  sf_warning("it=%d;",it+1);
#if defined _OPENMP && _DEBUG
    tic=omp_get_wtime();
#endif

    step_forward(u0,u1,vel,rho,fdcoef_d2,fdcoef_d1,nop,nzpad,nxpad,nypad);

    if (adj) { /* backward inject source wavelet */
      if (expl) {
        sf_seek(file_wav,(off_t)(nt-it-1)*sizeof(float),SEEK_SET);
        sf_floatread(ws,1,file_wav);
        if (sinc) sinc3d_inject1_with_vv(u0,ws[0],cssinc,vel);
        else      lint3d_inject1_with_vv(u0,ws[0],cslint,vel);
      } else {
        sf_seek(file_wav,(off_t)(nt-it-1)*ns*sizeof(float),SEEK_SET);
        sf_floatread(ws,ns,file_wav);
        if (sinc) sinc3d_inject_with_vv(u0,ws,cssinc,vel);
        else      lint3d_inject_with_vv(u0,ws,cslint,vel);
      }
    } else { /* forward inject source wavelet */
      if (expl) {
        sf_floatread(ws,1,file_wav);
        if (sinc) sinc3d_inject1_with_vv(u0,ws[0],cssinc,vel);
        else      lint3d_inject1_with_vv(u0,ws[0],cslint,vel);
      } else {
        sf_floatread(ws,ns,file_wav);
        if (sinc) sinc3d_inject_with_vv(u0,ws,cssinc,vel);
        else      lint3d_inject_with_vv(u0,ws,cslint,vel);
      }
    }

    /* apply abc */
    if (dabc) {
      if (hybrid) apply_abc(u0,u1,nz,nx,ny,nbd,abc,nop,damp);
      else {
        abcone3d_apply(u0,u1,nop,abc,fdm);
        sponge3d_apply(u0,spo,fdm);
        sponge3d_apply(u1,spo,fdm);
      }
    }

    /* loop over pointers */
    ptr_tmp = u0;  u0 = u1;  u1 = ptr_tmp;

    /* extract snapshot */
    if (snap && it%jsnap==0) {
      int fy = (floor)((sf_o(acy)-fdm->oypad)/fdm->dy);
      int jy = floor(sf_d(acy)/fdm->dy);
      float **ptr_slice;
      for (iy=0; iy<sf_n(acy); iy++) {
        ptr_slice = u0[fy+iy*jy];
        cut3d_slice(ptr_slice,oslice,fdm,acz,acx);
        sf_floatwrite(oslice[0],sf_n(acz)*sf_n(acx),file_wfl);
      }
    }

    /* extract receiver data */
    if (sinc) sinc3d_extract(u0,u_dat,crsinc);
    else      lint3d_extract(u0,u_dat,crlint);

    sf_floatwrite(u_dat,nr,file_dat);

#if defined _OPENMP && _DEBUG
    toc=omp_get_wtime();
    fprintf(stderr,"%5.2gs",(float)(toc-tic));
#endif
  }
#ifdef _OPENMP
  wall_clock_time_e = omp_get_wtime();
#else
  wall_clock_time_e = (double) clock() / CLOCKS_PER_SEC;
#endif
  if (verb)
    fprintf(stderr,"\nElapsed time: %lf s\n",wall_clock_time_e-wall_clock_time_s);

  free(**u0); free(*u0); free(u0);
  free(**u1); free(*u1); free(u1);
  free(**vel); free(*vel); free(vel);
  free(u_dat);
  free(ws);
  free(fdcoef_d2); free(fdcoef_d1);
  if (snap) { free(*oslice); free(oslice); }
  if(!cden) { free(**rho); free(*rho); free(rho); }
  if (hybrid) free(damp);
  free(src3d); free(rec3d);

  return 0;
}

