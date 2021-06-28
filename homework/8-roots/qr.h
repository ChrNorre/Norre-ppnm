void QR_decomp(gsl_matrix* A, gsl_matrix* R);
void QR_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void QR_backsub(gsl_matrix* R, gsl_vector* x);
void QR_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B, gsl_vector* x);
