#include<math.h>
#include"montecarlo.h"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#define get gsl_vector_get
#define set gsl_vector_set

void randomx(gsl_vector* a, gsl_vector* b, gsl_vector* x){
	for (int i=0; i<a->size; i++){
		double r = 1.0*rand()/RAND_MAX;
		set(x,i, get(a,i) +  r*(get(b,i) - get(a,i)));
	}
}

void plainmc(gsl_vector* a, gsl_vector* b, double (*f)(gsl_vector* x), int N, double* result, double* error){
	double V=1;
	for (int i=0; i<a->size; i++){
		V *= get(b,i)-get(a,i);
	}
	double sum1=0, sum2=0;
	double fx;
	gsl_vector* x = gsl_vector_alloc(a->size);
	for (int k=0; k<N; k++){
		randomx(a,b,x);
		fx = f(x);
		sum1 += fx;
		sum2 += fx*fx;
	}
	double avr = sum1/N;
	double var = sum2/N - avr*avr;
	*result = avr*V;
	*error = sqrt(var/N)*V;
	gsl_vector_free(x);
}


double corput(int n, int base){
	//calculate the n'th corput number in base 'base'
	double q = 0, bk = 1.0/base;
	while (n>0){ 
		q += (n % base)*bk; 
		n /= base;
		bk /= base;
	}
	return q;
}

void halton1(int n, gsl_vector* x){
	//calculates the n'th multidimensional halton number
	int base[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
	for(int i=0; i<x->size; i++){
		set(x,i,corput(n,base[i]));
	}
}

void halton2(int n, gsl_vector* x){
	//calculates the n'th multidimensional halton number
	//shifted base should function as another sequence for error calculation
	int base[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
	for(int i=0; i<x->size; i++){
		set(x,i,corput(n,base[i]));
	}
}

void qrandomx(gsl_vector* a, gsl_vector* b, gsl_vector* x){
	//scaling from quasi-random unit(ish)-vector to correct interval
	for (int i=0; i<a->size; i++){
		set(x,i, get(a,i) +  get(x,i)*(get(b,i) - get(a,i)));
	}

}


void quasimc(gsl_vector* a, gsl_vector* b, double (*f)(gsl_vector* x), int N, double* result, double* error){
	double V=1;
	for (int i=0; i<a->size; i++){
		V *= get(b,i)-get(a,i);
	}
	double sum1=0, sum2=0;
	gsl_vector* x1 = gsl_vector_alloc(a->size);
	gsl_vector* x2 = gsl_vector_alloc(a->size);
	for (int k=0; k<N; k++){
		halton1(k,x1);
		qrandomx(a,b,x1);
		halton2(k,x2);
		qrandomx(a,b,x2);
		sum1 += f(x1);
		sum2 += f(x2);
	}
	double I1 = V*sum1/N;
	double I2 = V*sum2/N;
	*result = (I1+I2)/2;
	*error = fabs(I2-I1);
	gsl_vector_free(x1);
	gsl_vector_free(x2);
}



