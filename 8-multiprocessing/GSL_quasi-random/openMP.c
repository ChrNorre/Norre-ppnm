#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<omp.h>
#include<gsl/gsl_qrng.h>

struct param {int N,res; gsl_qrng* q;};

void* pi(void* args);

int main(){
	
	int N = 3e7;
	int n1=N/3, n2=N/3, n3=N/3;

	gsl_qrng* qrng1 = gsl_qrng_alloc(gsl_qrng_halton, 2);
	gsl_qrng* qrng2 = gsl_qrng_alloc(gsl_qrng_halton, 2);
	gsl_qrng* qrng3 = gsl_qrng_alloc(gsl_qrng_halton, 2);

	struct param p1 = {.N=n1,.res=0,.q=qrng1};
	struct param p2 = {.N=n2,.res=0,.q=qrng2};
	struct param p3 = {.N=n3,.res=0,.q=qrng3};


#pragma omp parallel sections
	{
	#pragma omp section
		{
			pi((void*)&p1);
		}
	#pragma omp section
		{
			pi((void*)&p2);
		}
	#pragma omp section	
		{
			pi((void*)&p3);
		}
	}

	gsl_qrng_free(qrng1);
	gsl_qrng_free(qrng2);
	gsl_qrng_free(qrng3);
	

	double pi_est = (4.0*(p1.res + p2.res + p3.res))/N;
	printf("for N=%i, pi is estimated to be: %.12g\n",N,pi_est);

return 0;
}