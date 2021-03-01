#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<omp.h>
#include<gsl/gsl_qrng.h>

struct param {int N,res; gsl_qrng* q;};

void* pi(void* args);


int main(int argc, char** argv){

	int N, n1, n2, n3;
	double pi_est;
	struct param p1, p2, p3;

	gsl_qrng* qrng1 = gsl_qrng_alloc(gsl_qrng_halton, 2);
	gsl_qrng* qrng2 = gsl_qrng_alloc(gsl_qrng_halton, 2);
	gsl_qrng* qrng3 = gsl_qrng_alloc(gsl_qrng_halton, 2);

	if(argc>1) {
		for (int i= 1; i<argc; i++) {
	
			N = 3*pow(10,atof(argv[i]));
			n1=N/3;
			n2=N/3;
			n3=N/3;
			

			p1.N=n1; p1.res=0; p1.q=qrng1;
			p2.N=n2; p2.res=0; p2.q=qrng2;
			p3.N=n3; p3.res=0; p3.q=qrng3;


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
			

			pi_est = (4.0*(p1.res + p2.res + p3.res))/N;
			printf("%.16g %.16g\n",log10(N), log10(fabs(pi_est-M_PI)));
		}
	}

	gsl_qrng_free(qrng1);
	gsl_qrng_free(qrng2);
	gsl_qrng_free(qrng3);

return 0;
}