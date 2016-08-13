#ifndef EWEFD3D_KERNELS_H_TFI0CHLV
#define EWEFD3D_KERNELS_H_TFI0CHLV

void expand_cpu(float *a, float *b, int nb, int x_a, int x_b, int y_a, int y_b, int z_a, int z_b);

__global__ void computeRo(float *d_ro, float dt, int nxpad, int nzpad, int nyinterior);

__global__ void dispToStrain(int nxpad, int nylocal, int nzpad, float *d_uox, float *d_uoy, float *d_uoz, float *d_txx, float *d_tyy, float *d_tzz, float *d_txy, float *d_tyz, float *d_tzx, float idx, float idy, float idz);

__global__ void strainToStress(int gpuID, int nxpad, int nzpad, int nyinterior, float *d_c11, float *d_c12, float *d_c13, float *d_c22, float *d_c23, float *d_c33, float *d_c44, float *d_c55, float *d_c66, float *d_txx, float *d_tyy, float *d_tzz, float *d_txy, float *d_tyz, float *d_tzx);

__global__ void freeSurf(int gpuID, int nxpad, int nyinterior, int nzpad, int nb, float *d_tzz, float *d_tyz, float *d_tzx);

__global__ void lint3d_bell_gpu(int gpuID, int it, int nc, int ns, int c, int nbell, int nxpad, int nyinterior, int nzpad, float *d_uu, float *d_bell, int *d_jx, int *d_jz, int *d_jy, float *d_ww, float *d_Sw000, float *d_Sw001, float *d_Sw010, float *d_Sw011, float *d_Sw100, float *d_Sw101, float *d_Sw110, float *d_Sw111);


__global__ void stressToAccel(int nxpad, int nzpad, int nylocal, float idx, float idy, float idz, float *d_txx, float *d_tyy, float *d_tzz, float *d_txy, float *d_tzx, float *d_tyz, float *d_uax, float *d_uay, float *d_uaz);

__global__ void stepTime(int gpuID, int nxpad, int nyinterior, int nzpad, float *d_ro, float *d_uox, float *d_umx, float *d_uax, float *d_upx, float *d_uoy, float *d_umy, float *d_uay, float *d_upy, float *d_uoz, float *d_umz, float *d_uaz, float *d_upz);


__global__ void abcone3d_apply_XY(int gpuID, int nxpad, int nyinterior, int nzpad, float *d_uo, float *d_um, float *d_bzl, float *d_bzh);

__global__ void abcone3d_apply_ZY(int gpuID, int nxpad, int nyinterior, int nzpad, float *d_uo, float *d_um, float *d_bxl, float *d_bxh);

__global__ void abcone3d_apply_XZ_low(int nxpad, int nzpad, float *d_uo, float *d_um, float *d_byl);

__global__ void abcone3d_apply_XZ_high(int nxpad, int nylocal, int nzpad, float *d_uo, float *d_um, float *d_byh);

__global__ void sponge3d_apply_XY(int gpuID, float *d_uu, int nxpad, int nyinterior, int nzpad, int nb);

__global__ void sponge3d_apply_ZY(int gpuID, float *d_uu, int nxpad, int nyinterior, int nzpad, int nb, int nx, float spo);


__global__ void sponge3d_apply_XZ_low(float *d_uu, int nxpad, int nzpad, int nb, int nylocal);

__global__ void sponge3d_apply_XZ_high(float *d_uu, int nxpad, int nylocal, int nzpad, int nb);

__global__ void lint3d_extract_gpu(int gpuID, float *d_dd, int nr, int nxpad, int nyinterior, int nzpad, float *d_uoz, float *d_uox, float *d_uoy, int *d_Rjz, int *d_Rjx, int *d_Rjy, float *d_Rw000, float *d_Rw001, float *d_Rw010, float *d_Rw011, float *d_Rw100, float *d_Rw101, float *d_Rw110, float *d_Rw111);

__global__ void extract_gpu(int gpuID, float *d_dd, int nr, int nxpad, int nyinterior, int nzpad, float *d_uoz, float *d_uox, float *d_uoy, int *d_Rjz, int *d_Rjx, int *d_Rjy);

#endif /* end of include guard: EWEFD3D_KERNELS_H_TFI0CHLV */
