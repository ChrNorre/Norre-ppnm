#include<gsl/gsl_vector.h>
typedef struct { int n; double(*f)(double); gsl_vector* params; } ann;
ann*   ann_alloc   (int n,double(*f)(double));
void   ann_free    (ann* network);
double ann_response(ann* network,double x);
void   ann_train   (ann* network,gsl_vector* xs,gsl_vector* ys);
double ann_integ(ann* network, double (*F)(double), double z);
double ann_deriv(ann* network, double (*df)(double), double z);