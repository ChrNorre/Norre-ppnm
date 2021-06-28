#include<gsl/gsl_vector.h>
void randomx(gsl_vector* a, gsl_vector* b, gsl_vector* x);
void plainmc(gsl_vector* a, gsl_vector* b, double (*f)(gsl_vector* x), int N, double* result, double* error);
void quasimc(gsl_vector* a, gsl_vector* b, double (*f)(gsl_vector* x), int N, double* result, double* error);
