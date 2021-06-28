#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include"qr.h"
#include"rootfinder.h"

void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps){
	int n = x->size;
	gsl_matrix* J = gsl_matrix_alloc(n,n);	//Jacobian
	gsl_matrix* R = gsl_matrix_alloc(n,n);	//extra matrix for solving linear system
	gsl_vector* fx = gsl_vector_alloc(n);	//function value
	gsl_vector* fx_lambda = gsl_vector_alloc(n);	//function value at x+lambda*dx
	gsl_vector* dx = gsl_vector_alloc(n);	//step size
	gsl_vector* x_lambda = gsl_vector_alloc(n);	//x value at x+lambda*dx
	double delx = 1e-6;
	double lambda;


	//iterate		A do-while whould have been cleaner
	while(1){

		f(x,fx); 	//calculates f(x) and saves in fx

		//calculating Jacobian
		for (int j = 0; j < n; j++){
			gsl_vector_set(x, j, gsl_vector_get(x,j)+delx); //x -> x+delx
			f(x,fx_lambda);
			gsl_blas_daxpy(-1, fx, fx_lambda); //fx_lambda = f(x+delx)-f(x) = df

			for (int k = 0; k < n; k++){
				gsl_matrix_set(J,k,j, 1.0*gsl_vector_get(fx_lambda, k)/delx);	//J(k,j) = df(k)/delx
			}
			gsl_vector_set(x, j, gsl_vector_get(x,j)-delx);
		}


		//Solving jacobian to get step size dx
		gsl_vector_scale(fx, -1);
		QR_decomp(J, R);
		QR_solve(J, R, fx, dx);
		gsl_vector_scale(fx, -1);


		lambda = 1.0;
		while(1){
			//scaling stepsize by lambda to avoid overshooting the root

			gsl_vector_memcpy(x_lambda, x);		//x -> x_lambda
			gsl_blas_daxpy(lambda, dx, x_lambda);	//x_lambda = x + lambda*dx
			f(x_lambda,fx_lambda);				//calculates f(x_lambda)

			if ((gsl_blas_dnrm2(fx_lambda) < (1-lambda/2)*gsl_blas_dnrm2(fx)) || (lambda < 1.0/64)){
				// if our current step gives a significant drop in the size of the function
				// OR if lambda gets too small, it signals problematic point, step anyways
				break;
			}
			lambda *= 0.5;
		}
		gsl_vector_memcpy(x, x_lambda);		//apply step size, x_lambda -> x

		//Until our step dx becomes too small
		//or our function at x becomes sufficiently small
		if((gsl_blas_dnrm2(dx)<delx) || (gsl_blas_dnrm2(fx)<eps)){
			break;
		}
	}

	//when out of while, root is saved in x
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(fx);
	gsl_vector_free(fx_lambda);
	gsl_vector_free(dx);
	gsl_vector_free(x_lambda);
}