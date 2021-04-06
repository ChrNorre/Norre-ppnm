void least_squares(gsl_vector* xs, gsl_vector* ys, gsl_vector* dys, double (*f)(int, double), gsl_vector* cs, gsl_matrix* Sigma);
double eval_fun(double (*f)(int, double), gsl_vector* cs, double x);