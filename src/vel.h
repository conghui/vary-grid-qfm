#ifndef VEL_H_BLY3XMED
#define VEL_H_BLY3XMED


typedef struct {
  int   nz, nx, ny;
  float oz, ox, oy;
  float dz, dx, dy;
  float ***dat;
  float ***vgamma;
} vel_t;

#endif /* end of include guard: VEL_H_BLY3XMED */
