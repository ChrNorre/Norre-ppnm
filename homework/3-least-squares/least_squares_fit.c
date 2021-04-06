#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include"qr.h"


void least_squares(gsl_vector* xs, gsl_vector* ys, gsl_vector* dys, double (*f)(int, double), gsl_vector* cs, gsl_matrix* Sigma){
	//Linear Fitting
	//takes vector of x-values, y-values, error in y-values, function to fit to (f)
	//vector of fitting coefficients, cs, (initial value irrelevant, length is the most important)
	//Optional covariance matrix Sigma (square of same size as cs) (initial value irrelevant), Replace with NULL, if covariance is unwanted

	int n = xs->size;
	int k = cs->size;

	gsl_matrix* A = gsl_matrix_alloc(n,k);
	gsl_matrix* R = gsl_matrix_alloc(k,k);
	gsl_vector* b = gsl_vector_alloc(n);

	double res;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < k; j++){
			res = f(j,gsl_vector_get(xs,i));
			gsl_matrix_set(A,i,j, res/gsl_vector_get(dys,i));
		}
		gsl_vector_set(b,i, gsl_vector_get(ys,i)/gsl_vector_get(dys,i));
	}

	QR_decomp(A, R);

	QR_solve(A, R, b, cs);

	//covariance matrix
	if (Sigma != NULL){
		//Initializa as identity matrix
		for (int i = 0; i < k; i++){
			for (int j = 0; j < k; j++){
				if (i==j){
					gsl_matrix_set(Sigma, i, j, 1);
				}
				else {
					gsl_matrix_set(Sigma, i, j, 0);
				}
			}
		}

		//calculates R^-1 and saves in Sigma
		gsl_vector_view x;
		for (int i=0; i<k; i++){
			x = gsl_matrix_column(Sigma,i);
			QR_backsub(R, &x.vector);
		}

		gsl_matrix_memcpy(R, Sigma);

		gsl_blas_dtrmm(CblasRight, CblasUpper, CblasTrans, CblasNonUnit, 1, R, Sigma); //Sigma = 1*Sigma*(R^T)
		//unsure about CblasNonUnit		
	}

	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_vector_free(b);
}


double eval_fun(double (*f)(int, double), gsl_vector* cs, double x){
	int k = cs->size;
	double res = 0;
	for (int i = 0; i < k; i++){
		res += gsl_vector_get(cs,i) * f(i,x);
	}
	return res;
}