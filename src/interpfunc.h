#ifndef INTERPFUNC_H_CJ7TXNGP
#define INTERPFUNC_H_CJ7TXNGP

void make_sinc_table(float *sinc_table, int ntab, int lsinc);

void sinc_interp3d_1(const float *input, float *output, const float *sinc_table,
    int ntab, int lsinc,
    int n1, float o1, float d1,
    int n2, float o2, float d2,
    int n3, float o3, float d3,
    int n1out, float o1out, float d1out,
    int n2out, float o2out, float d2out,
    int n3out, float o3out, float d3out);

#endif /* end of include guard: INTERPFUNC_H_CJ7TXNGP */
