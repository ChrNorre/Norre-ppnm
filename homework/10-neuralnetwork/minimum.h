void num_gradient(double F(gsl_vector* x), gsl_vector* x, gsl_vector* gx);
void qnewton(
	double F(gsl_vector* x), /* objective function */
	gsl_vector* x, /* on input: starting point, on exit: approximation to root */
	double eps, /* accuracy goal, on exit |gradient| should be <eps */
	int* counter
);