static void sinc_interp3d_1_(const fdm3d &oldfdm, const fdm3d &newfdm, const float *sinc_table, int ntab, int lsinc, float ***u, float ***nu)
{
  sinc_interp3d_1(u[0][0], nu[0][0], sinc_table, ntab, lsinc,
    oldfdm->nzpad, oldfdm->oz, oldfdm->dz,  /* old */
    oldfdm->nxpad, oldfdm->ox, oldfdm->dx,
    oldfdm->nypad, oldfdm->oy, oldfdm->dy,
    newfdm->nzpad, newfdm->oz, newfdm->dz,  /* new */
    newfdm->nxpad, newfdm->ox, newfdm->dx,
    newfdm->nypad, newfdm->oy, newfdm->dy);
}

static void interp_host_umo_patch_2(const fdm3d &oldfdm, const fdm3d &newfdm, const float *sinc_table, int ntab, int lsinc, float ***&h_umx, float ***&h_uox,  float ***&h_umy,  float ***&h_uoy,  float ***&h_umz,  float ***&h_uoz)
{
  struct timeval start, stop;
  gettimeofday(&start, NULL);

  float ***n_umx = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  float ***n_umy = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  float ***n_umz = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  float ***n_uox = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  float ***n_uoy = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  float ***n_uoz = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);

  sinc_interp3d_1_(oldfdm, newfdm, sinc_table, ntab, lsinc, h_umx, n_umx);
  sinc_interp3d_1_(oldfdm, newfdm, sinc_table, ntab, lsinc, h_umy, n_umy);
  sinc_interp3d_1_(oldfdm, newfdm, sinc_table, ntab, lsinc, h_umz, n_umz);
  sinc_interp3d_1_(oldfdm, newfdm, sinc_table, ntab, lsinc, h_uox, n_uox);
  sinc_interp3d_1_(oldfdm, newfdm, sinc_table, ntab, lsinc, h_uoy, n_uoy);
  sinc_interp3d_1_(oldfdm, newfdm, sinc_table, ntab, lsinc, h_uoz, n_uoz);

  free(**h_umx); free(*h_umx); free(h_umx); h_umx = n_umx;
  free(**h_umy); free(*h_umy); free(h_umy); h_umy = n_umy;
  free(**h_umz); free(*h_umz); free(h_umz); h_umz = n_umz;
  free(**h_uox); free(*h_uox); free(h_uox); h_uox = n_uox;
  free(**h_uoy); free(*h_uoy); free(h_uoy); h_uoy = n_uoy;
  free(**h_uoz); free(*h_uoz); free(h_uoz); h_uoz = n_uoz;

  gettimeofday(&stop, NULL);
  float elapse = stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) * 1e-6;
  fprintf(stderr, "time for interpolate wavefields: %.2f\n", elapse);
}

static void interp_host_den_vel_patch_2(const fdm3d &oldfdm, const fdm3d &newfdm, const fdm3d &fullfdm, const float *sinc_table, int ntab, int lsinc, float ***full_h_ro, float ***full_h_c11, float ***full_h_c22, float ***full_h_c33, float ***full_h_c44, float ***full_h_c55, float ***full_h_c66, float ***full_h_c12, float ***full_h_c13, float ***full_h_c23, float ***&h_ro, float ***&h_c11, float ***&h_c22, float ***&h_c33, float ***&h_c44, float ***&h_c55, float ***&h_c66, float ***&h_c12, float ***&h_c13, float ***&h_c23)
{
  if (checksame(oldfdm, newfdm)) {
    sf_warning("old and new dimensions for den/vel are the same, don't need interpolatin");
    return;
  }

  struct timeval start, stop;
  gettimeofday(&start, NULL);

  free(**h_ro); free(*h_ro); free(h_ro); h_ro = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c11); free(*h_c11); free(h_c11); h_c11 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c22); free(*h_c22); free(h_c22); h_c22 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c33); free(*h_c33); free(h_c33); h_c33 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c44); free(*h_c44); free(h_c44); h_c44 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c55); free(*h_c55); free(h_c55); h_c55 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c66); free(*h_c66); free(h_c66); h_c66 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c12); free(*h_c12); free(h_c12); h_c12 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c13); free(*h_c13); free(h_c13); h_c13 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c23); free(*h_c23); free(h_c23); h_c23 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);


  sinc_interp3d_1_(fullfdm, newfdm, sinc_table, ntab, lsinc, full_h_ro, h_ro);
  sinc_interp3d_1_(fullfdm, newfdm, sinc_table, ntab, lsinc, full_h_c11, h_c11);
  sinc_interp3d_1_(fullfdm, newfdm, sinc_table, ntab, lsinc, full_h_c22, h_c22);
  sinc_interp3d_1_(fullfdm, newfdm, sinc_table, ntab, lsinc, full_h_c33, h_c33);
  sinc_interp3d_1_(fullfdm, newfdm, sinc_table, ntab, lsinc, full_h_c44, h_c44);
  sinc_interp3d_1_(fullfdm, newfdm, sinc_table, ntab, lsinc, full_h_c55, h_c55);
  sinc_interp3d_1_(fullfdm, newfdm, sinc_table, ntab, lsinc, full_h_c66, h_c66);
  sinc_interp3d_1_(fullfdm, newfdm, sinc_table, ntab, lsinc, full_h_c12, h_c12);
  sinc_interp3d_1_(fullfdm, newfdm, sinc_table, ntab, lsinc, full_h_c13, h_c13);
  sinc_interp3d_1_(fullfdm, newfdm, sinc_table, ntab, lsinc, full_h_c23, h_c23);

  gettimeofday(&stop, NULL);
  float elapse = stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) * 1e-6;
  fprintf(stderr, "time for interpolate density and velocities: %.2f\n", elapse);
}

static void interp_host_den_vel_patch(const fdm3d &oldfdm, const fdm3d &newfdm, const fdm3d &fullfdm, float ***full_h_ro, float ***full_h_c11, float ***full_h_c22, float ***full_h_c33, float ***full_h_c44, float ***full_h_c55, float ***full_h_c66, float ***full_h_c12, float ***full_h_c13, float ***full_h_c23, float ***&h_ro, float ***&h_c11, float ***&h_c22, float ***&h_c33, float ***&h_c44, float ***&h_c55, float ***&h_c66, float ***&h_c12, float ***&h_c13, float ***&h_c23)
{
  if (checksame(oldfdm, newfdm)) {
    sf_warning("old and new dimensions for den/vel are the same, don't need interpolatin");
    return;
  }

  struct timeval start, stop;
  gettimeofday(&start, NULL);

  free(**h_ro); free(*h_ro); free(h_ro); h_ro = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c11); free(*h_c11); free(h_c11); h_c11 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c22); free(*h_c22); free(h_c22); h_c22 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c33); free(*h_c33); free(h_c33); h_c33 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c44); free(*h_c44); free(h_c44); h_c44 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c55); free(*h_c55); free(h_c55); h_c55 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c66); free(*h_c66); free(h_c66); h_c66 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c12); free(*h_c12); free(h_c12); h_c12 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c13); free(*h_c13); free(h_c13); h_c13 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  free(**h_c23); free(*h_c23); free(h_c23); h_c23 = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);


  interp_den_vel_(
      full_h_ro, full_h_c11, full_h_c22, full_h_c33,
      full_h_c44, full_h_c55, full_h_c66, full_h_c12,
      full_h_c13, full_h_c23,
      h_ro, h_c11, h_c22, h_c33, h_c44,
      h_c55, h_c66, h_c12, h_c13, h_c23,
    fullfdm->nzpad, fullfdm->oz, fullfdm->dz,  /* old */
    fullfdm->nxpad, fullfdm->ox, fullfdm->dx,
    fullfdm->nypad, fullfdm->oy, fullfdm->dy,
    newfdm->nzpad, newfdm->oz, newfdm->dz,  /* new */
    newfdm->nxpad, newfdm->ox, newfdm->dx,
    newfdm->nypad, newfdm->oy, newfdm->dy);

  gettimeofday(&stop, NULL);
  float elapse = stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) * 1e-6;
  fprintf(stderr, "time for interpolate density and velocities: %.2f\n", elapse);
}

static void interp_host_umo_patch(const fdm3d &oldfdm, const fdm3d &newfdm, float ***&h_umx, float ***&h_uox,  float ***&h_umy,  float ***&h_uoy,  float ***&h_umz,  float ***&h_uoz)
{
  struct timeval start, stop;
  gettimeofday(&start, NULL);

  float ***n_umx = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  float ***n_umy = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  float ***n_umz = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  float ***n_uox = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  float ***n_uoy = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  float ***n_uoz = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);

  interp_wavefield_(
    h_umx, h_uox,  h_umy,  h_uoy,  h_umz,  h_uoz,
    n_umx, n_uox,  n_umy,  n_uoy,  n_umz,  n_uoz,
    oldfdm->nzpad, oldfdm->oz, oldfdm->dz,  /* old */
    oldfdm->nxpad, oldfdm->ox, oldfdm->dx,
    oldfdm->nypad, oldfdm->oy, oldfdm->dy,
    newfdm->nzpad, newfdm->oz, newfdm->dz,  /* new */
    newfdm->nxpad, newfdm->ox, newfdm->dx,
    newfdm->nypad, newfdm->oy, newfdm->dy);

  free(**h_umx); free(*h_umx); free(h_umx); h_umx = n_umx;
  free(**h_umy); free(*h_umy); free(h_umy); h_umy = n_umy;
  free(**h_umz); free(*h_umz); free(h_umz); h_umz = n_umz;
  free(**h_uox); free(*h_uox); free(h_uox); h_uox = n_uox;
  free(**h_uoy); free(*h_uoy); free(h_uoy); h_uoy = n_uoy;
  free(**h_uoz); free(*h_uoz); free(h_uoz); h_uoz = n_uoz;

  gettimeofday(&stop, NULL);
  float elapse = stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) * 1e-6;
  fprintf(stderr, "time for interpolate wavefields: %.2f\n", elapse);
}

static void interp_host_umo(const fdm3d &oldfdm, const fdm3d &newfdm, float ***&h_umx, float ***&h_uox,  float ***&h_umy,  float ***&h_uoy,  float ***&h_umz,  float ***&h_uoz)
{
  interp_wavefield(oldfdm, newfdm, h_umx);
  interp_wavefield(oldfdm, newfdm, h_umy);
  interp_wavefield(oldfdm, newfdm, h_umz);
  interp_wavefield(oldfdm, newfdm, h_uox);
  interp_wavefield(oldfdm, newfdm, h_uoy);
  interp_wavefield(oldfdm, newfdm, h_uoz);
}

static void interp_host_den_vel(const fdm3d &oldfdm, const fdm3d &newfdm, float ***full_h_ro, float ***full_h_c11, float ***full_h_c22, float ***full_h_c33, float ***full_h_c44, float ***full_h_c55, float ***full_h_c66, float ***full_h_c12, float ***full_h_c13, float ***full_h_c23, float ***&h_ro, float ***&h_c11, float ***&h_c22, float ***&h_c33, float ***&h_c44, float ***&h_c55, float ***&h_c66, float ***&h_c12, float ***&h_c13, float ***&h_c23)
{
  struct timeval start, stop;
  gettimeofday(&start, NULL);

  interp_den_vel(oldfdm, newfdm, full_h_ro, h_ro);
  interp_den_vel(oldfdm, newfdm, full_h_c11, h_c11);
  interp_den_vel(oldfdm, newfdm, full_h_c22, h_c22);
  interp_den_vel(oldfdm, newfdm, full_h_c33, h_c33);
  interp_den_vel(oldfdm, newfdm, full_h_c44, h_c44);
  interp_den_vel(oldfdm, newfdm, full_h_c55, h_c55);
  interp_den_vel(oldfdm, newfdm, full_h_c66, h_c66);
  interp_den_vel(oldfdm, newfdm, full_h_c12, h_c12);
  interp_den_vel(oldfdm, newfdm, full_h_c13, h_c13);
  interp_den_vel(oldfdm, newfdm, full_h_c23, h_c23);

  gettimeofday(&stop, NULL);
  float elapse = stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) * 1e-6;
  fprintf(stderr, "time for interpolate density and velocities: %.2f\n", elapse);
}

static void interp_wavefield(const fdm3d &oldfdm, const fdm3d &newfdm, float ***&u)
{
  if (checksame(oldfdm, newfdm)) {
    sf_warning("old and new dimensions for wavefields are the same, don't need interpolatin");
    return;
  }

  float ***nu = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  interpfield_(u, nu, false,
    oldfdm->nzpad, oldfdm->oz, oldfdm->dz,  /* old */
    oldfdm->nxpad, oldfdm->ox, oldfdm->dx,
    oldfdm->nypad, oldfdm->oy, oldfdm->dy,
    newfdm->nzpad, newfdm->oz, newfdm->dz,  /* new */
    newfdm->nxpad, newfdm->ox, newfdm->dx,
    newfdm->nypad, newfdm->oy, newfdm->dy);

  free(**u); free(*u); free(u);
  u = nu;
}


static void interp_den_vel(const fdm3d &oldfdm, const fdm3d &newfdm, float ***oldf, float ***&newf)
{
  if (checksame(oldfdm, newfdm)) {
    sf_warning("old and new dimensions for den/vel are the same, don't need interpolatin");
    return;
  }
  float ***nu = sf_floatalloc3(newfdm->nzpad, newfdm->nxpad, newfdm->nypad);
  interpfield_(oldf, nu, true,
    oldfdm->nzpad, oldfdm->oz, oldfdm->dz,  /* old */
    oldfdm->nxpad, oldfdm->ox, oldfdm->dx,
    oldfdm->nypad, oldfdm->oy, oldfdm->dy,
    newfdm->nzpad, newfdm->oz, newfdm->dz,  /* new */
    newfdm->nxpad, newfdm->ox, newfdm->dx,
    newfdm->nypad, newfdm->oy, newfdm->dy);

  free(**newf); free(*newf); free(newf);
  newf = nu;
}
static void set_output_wfd(sf_file &Fwfl, sf_axis &at, const sf_axis &az, const sf_axis &ax, const sf_axis &ay, const sf_axis &ac, int nt, float dt, int jsnap, bool verb)
{
  int ntsnap=0;
  for(int it=0; it<nt; it++) {
    if(it%jsnap==0 && it != 0) ntsnap++;
  }
  sf_setn(at,  ntsnap);
  sf_setd(at,dt*jsnap);
  if(verb) sf_raxa(at);

  sf_oaxa(Fwfl,az,1);
  sf_oaxa(Fwfl,ax,2);
  sf_oaxa(Fwfl,ay,3);
  sf_oaxa(Fwfl,ac, 4);
  sf_oaxa(Fwfl,at, 5);
}
