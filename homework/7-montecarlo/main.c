#include"montecarlo.h"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#define get gsl_vector_get
#define set gsl_vector_set


int main(){

	//classic example of integrating a constant unit disc inside a square
	//should give pi
	double disc(gsl_vector* x){
		if(gsl_blas_dnrm2(x) < 1) {return 1;}
		else {return 0;}
	}
	int n = 2;

	gsl_vector* a1 = gsl_vector_alloc(n);
	gsl_vector* b1 = gsl_vector_alloc(n);
	for (int i = 0; i<n; i++){
		set(a1,i,-1);
		set(b1,i,1);
	}

	double result, error;
	int Npow = 6;
	plainmc(a1, b1, &disc, pow(10,Npow), &result, &error);
	printf("\nIntegral over unit disc for N=1e%d  =  %5lg +/- %5lg\n",Npow,result,error);
	gsl_vector_free(a1);
	gsl_vector_free(b1);




	//troublesome integral
	//limits are from 0 to PI in all dimensions, the more M_PI i use, the more precision i lose
	//I let limits vary from 0 to 3.15, so volume is more exact
	//auxilliary function is just checking if each dimension is less than M_PI
	double gamma_int(gsl_vector* x){
		double xp = get(x,0);
		double yp = get(x,1);
		double zp = get(x,2);
		if(((xp<M_PI) && (yp<M_PI)) && (zp<M_PI)) {
			double r = pow(M_PI,3)*(1 - cos(xp)*cos(yp)*cos(zp));
			return 1.0/r;
		}
		return 0;		
	}
	n = 3;		//3 dim, (x,y,z)
	gsl_vector* a2 = gsl_vector_alloc(n);
	gsl_vector* b2 = gsl_vector_alloc(n);
	for (int i = 0; i<n; i++){
		set(a2,i,0);
		set(b2,i,3.15);
	}

	Npow = 6;
	plainmc(a2, b2, &gamma_int, pow(10,Npow), &result, &error);
	printf("\nIntegral over complicated function for N=1e%d  =  %.6lg +/- %.3lg\t |diff|=%.3lg +/- %.3lg\n",Npow,result,error,fabs(1.39320392968567685-result),error);
	gsl_vector_free(a2);
	gsl_vector_free(b2);




	//f(r) = |x|*exp(-r*r) over unit sphere in 8D, Auxilliary is to go from box to sphere
	// I say this is an appropriate use of monte carlo, 
	//		since no one wants to do 8D spherical coordinate transformations
	double test3(gsl_vector* x){
		double r = gsl_blas_dnrm2(x);
		if(r<1) {
			return fabs(get(x,0))*exp(-r*r);
		}
		return 0;		
	}
	n = 8;		//8 dim
	gsl_vector* a3 = gsl_vector_alloc(n);
	gsl_vector* b3 = gsl_vector_alloc(n);
	for (int i = 0; i<n; i++){
		set(a3,i,-1);
		set(b3,i,+1);
	}

	Npow = 6;
	plainmc(a3, b3, &test3, pow(10,Npow), &result, &error);
	printf("\nIntegral of |x|*exp(-r*r) over 8D unit sphere for N=1e%d  =  %.6lg +/- %.3lg\n",Npow,result,error);
	gsl_vector_free(a3);
	gsl_vector_free(b3);




	printf("\nComparing scaling of error\n See produced figure\n");
	//using the triple cosine integral

	n = 3;		//3 dim, (x,y,z)
	gsl_vector* a5 = gsl_vector_alloc(n);
	gsl_vector* b5 = gsl_vector_alloc(n);
	for (int i = 0; i<n; i++){
		set(a5,i,-1);
		set(b5,i,+1);	
	}

	double resp, errp, resq, errq;
	
	FILE* data = fopen("out.data.txt","w");
	for (int i=1; i<7; i++){
		plainmc(a5, b5, &test3, pow(10,i), &resp, &errp);
		quasimc(a5, b5, &test3, pow(10,i), &resq, &errq);
		fprintf(data,"%d %lg %lg\n",i, log(errp), log(errq));
	}		
	fclose(data);


	gsl_vector_free(a5);
	gsl_vector_free(b5);




	return 0;
}