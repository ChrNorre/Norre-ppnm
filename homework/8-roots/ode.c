#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include"ode.h"

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
){
	f(t, yt, k0); //gives k0 (slope at y(t))
	
	gsl_vector_memcpy(yh,yt);	
	gsl_blas_daxpy(0.5*h, k0, yh); // yh = 0.5*h*k0 + yt    half-step
	f(t+0.5*h,yh,k_half); //gives k_half (slope at half-step)

	gsl_vector_memcpy(yh,yt);
	f(t+h,yh,k1); //gives k1 (slope at full-step)

	gsl_blas_daxpy(4, k_half, k0);  //k0 = k0 + 4*k_half
	gsl_blas_daxpy(1, k1, k0);  //k0 = k0 + k1 =      k0 + 4 k_half + k1
	gsl_vector_scale(k0,1.0/6);		//k = k0     equation 16 in chapter "odes"

	gsl_blas_daxpy(h, k0, yh); // yh = h*k0 + yt    full-step by final slope

	//error should be h*(k-k_half) (third order method minus second order method)
	gsl_vector_memcpy(err,k0);
	gsl_blas_daxpy(-1,k_half,err);   // h*(k - k_half) = y - y(half-step slope)
	gsl_vector_scale(err,h);	// difference in result by using half-step slope vs third order is approximate error
}


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
){ /* calculate y(b) */

	int n = ya->size;
	gsl_vector* yt = gsl_vector_alloc(n);
	gsl_vector* yh = gsl_vector_alloc(n);
	gsl_vector* err = gsl_vector_alloc(n);
	gsl_vector* k0 = gsl_vector_alloc(n);
	gsl_vector* k_half = gsl_vector_alloc(n);
	gsl_vector* k1 = gsl_vector_alloc(n);
	
	double t = a;
	double hi = h;
	double ei, taui, h_accepted;
	gsl_vector_memcpy(yt, ya);

	if (path != NULL){
		fprintf(path,"%10g",t);
		for (int i = 0; i < n; i++){
			fprintf(path," %10g",gsl_vector_get(yt,i));
		}
	}

	while (t<b){ //while not done
		if (b-t < hi){ //if close to end
			hi = b-t;	//jump to end
		}
		do{	//do a step
			rkstep23(f, t, yt, hi, yh, err, k0, k_half, k1);

			ei = gsl_blas_dnrm2(err);
			taui = (eps*gsl_blas_dnrm2(yh) + acc) * sqrt(1.0*hi/(b-a));
			h_accepted = hi;
			hi = hi*pow(taui/ei,0.25)*0.95; 	//update step-size

		} while(ei>taui);	//if step was too big, repeat with new step but same t and yt
		
		gsl_vector_memcpy(yt, yh);	//if step was accepted, yh -> yt and update t, repeat until t=b
		t += h_accepted;

		if (path != NULL){
			fprintf(path,"\n%10g",t);
			for (int i = 0; i < n; i++){
				fprintf(path," %10g",gsl_vector_get(yt,i));
			}
		}
	}


	gsl_vector_memcpy(yb, yh);
	gsl_vector_free(yt);
	gsl_vector_free(yh);
	gsl_vector_free(err);
	gsl_vector_free(k0);
	gsl_vector_free(k_half);
	gsl_vector_free(k1);
}