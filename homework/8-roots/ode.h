void rkstep23(
	void (*f)(double, gsl_vector*, gsl_vector*), /* the f from dy/dt=f(t,y) */
	double t,              /* the current value of the variable */
	gsl_vector* yt,            /* the current value y(t) of the sought function */
	double h,              /* the step to be taken */
	gsl_vector* yh,             /* output: y(t+h) */
	gsl_vector* err,             /* output: error estimate */
	gsl_vector* k0,
	gsl_vector* k_half,			//k-vectors are just for calculations (input output irrelevant)
	gsl_vector* k1
);

void driver(
	void (*f)(double,gsl_vector*,gsl_vector*), /* right-hand-side of dy/dt=f(t,y) */
	double a,                     /* the start-point a */
	gsl_vector* ya,                     /* y(a) */
	double b,                     /* the end-point of the integration */
	gsl_vector* yb,                     /* y(b) to be calculated */
	double h,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps,                    /* relative accuracy goal */
	FILE* path                     /*FILE to write path */
);