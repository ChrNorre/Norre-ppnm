int binsearch_vec(gsl_vector* x, double z);
void find_ps(gsl_vector* x, gsl_vector* y, gsl_vector* p);
void cub_sub(gsl_vector* x, gsl_vector* y, gsl_vector* p, gsl_vector* c, gsl_vector* d);
double cub_sub_eval(gsl_vector* x, gsl_vector* y, gsl_vector* p, gsl_vector* c, gsl_vector* d, double z);
double cub_sub_deriv(gsl_vector* x, gsl_vector* y, gsl_vector* p, gsl_vector* c, gsl_vector* d, double z);
double cub_sub_integ(gsl_vector* x, gsl_vector* y, gsl_vector* p, gsl_vector* c, gsl_vector* d, double z);
