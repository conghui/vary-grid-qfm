#ifndef CHECK_H_NMIEIDTY
#define CHECK_H_NMIEIDTY

float fsum1(const float *a, int n);
int isum1(const int *a, int n);
void write3df(const char *file, float ***dat, int n1, int n2, int n3);
void write2di(const char *file, int **dat, int n1, int n2);
void write2df(const char *file, float **dat, int n1, int n2);
void write1di(const char *file, int *dat, int n1);


#endif /* end of include guard: CHECK_H_NMIEIDTY */
