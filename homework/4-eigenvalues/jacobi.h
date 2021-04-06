void timesJ(gsl_matrix* A, int p, int q, double theta);
void Jtimes(gsl_matrix* A, int p, int q, double theta);
void jacobi_diag(gsl_matrix* A, gsl_matrix* V);
int equal(double a, double b, double tau, double epsilon);