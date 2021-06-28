#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>

void num_gradient(double F(gsl_vector* x), gsl_vector* x, gsl_vector* gx){
	//calculates the numerical gradient of function F at point x and saves result in vector gx
	double delx = 1e-6;
	for (int i=0; i<gx->size; i++){
		double fx = F(x);
		gsl_vector_set(x,i,gsl_vector_get(x,i) + delx);
		double df = F(x) - fx;
		gsl_vector_set(gx,i,df/delx);
		gsl_vector_set(x,i,gsl_vector_get(x,i) - delx);
	}
}


void qnewton(
	double F(gsl_vector* x), /* objective function */
	gsl_vector* x, /* on input: starting point, on exit: approximation to root */
	double eps, /* accuracy goal, on exit |gradient| should be <eps */
	int* counter
){
	int n = x->size;
	gsl_matrix* B = gsl_matrix_alloc(n,n);
	gsl_matrix* delB = gsl_matrix_alloc(n,n);
	gsl_vector* dx = gsl_vector_calloc(n);	//step size
	gsl_vector* gx = gsl_vector_alloc(n);	//gradient
	gsl_vector* xs = gsl_vector_alloc(n);	//x + s = x + lambda*dx
	gsl_vector* u = gsl_vector_alloc(n);	//s = lambda*dx
	gsl_vector* y = gsl_vector_alloc(n);
	double delx = 1e-7;
	double lambda;
	double alpha = 1e-3;
	double res;
	*counter = 0;

	gsl_matrix_set_identity(B);

	
	//iterate	
	while(1){


		//initial step, equation (6)
		num_gradient(F,x,gx);
		gsl_blas_dgemv(CblasNoTrans, -1, B, gx, 0, dx); // dx = -B*gx


		//backtracking linesearch
		lambda = 1.0;
		while(1){
			//scaling stepsize by lambda to avoid overshooting
			gsl_vector_memcpy(xs, x);		//x -> s
			gsl_blas_daxpy(lambda, dx, xs);	//xs = x + lambda*dx = x + s

			gsl_blas_ddot(dx, gx, &res);   //dx^T * gx   ->   res

			if (F(xs) < F(x) + alpha*lambda*res){
				break;
			}
			if (lambda < 1.0/64){
				gsl_matrix_set_identity(B);
				break;
			}
			lambda *= 0.5;
		}

		gsl_vector_memcpy(u, dx);
		gsl_vector_scale(u,lambda);	// u = s = dx*lambda



		//updating inverse Hessian
		//delB = u uT / (uT y)

		num_gradient(F,xs,y);
		gsl_blas_daxpy(-1, gx, y);	// y = grad(F(x+s)) - grad(F(x))

		gsl_blas_dgemv(CblasNoTrans, -1, B, y, 1, u);	//u = s - By

		gsl_blas_ddot(u, y, &res);

		gsl_matrix_set_zero(delB);
		gsl_blas_dger(1/res, u, u, delB);     //delB = u uT / (uT y)

		if (abs(res) > 1e-6){
			gsl_matrix_add(B, delB);      //update  B -> B+delB
		}


		gsl_vector_memcpy(x, xs);		//apply step size, x+s -> new x

		//Until our step dx becomes too small
		//or our function at x becomes sufficiently small
		num_gradient(F,x,gx);
		*counter += 1;

		if((gsl_blas_dnrm2(dx)<delx) || (gsl_blas_dnrm2(gx)<eps)){
			break;
		}
	}

	//when out of while-loop, minimization saved in x

	gsl_matrix_free(B);
	gsl_matrix_free(delB);
	gsl_vector_free(dx);
	gsl_vector_free(gx);
	gsl_vector_free(xs);
	gsl_vector_free(u);
	gsl_vector_free(y);
}