#ifndef STEP_FORWARD_H_PL2PUYAU
#define STEP_FORWARD_H_PL2PUYAU

void step_forward(float*** restrict u0, float*** restrict u1,
    float*** restrict vel, float*** restrict rho,
    float* restrict fdcoef_d2, float* restrict fdcoef_d1,
    int nop, int nzpad, int nxpad, int nypad);


#endif /* end of include guard: STEP_FORWARD_H_PL2PUYAU */
