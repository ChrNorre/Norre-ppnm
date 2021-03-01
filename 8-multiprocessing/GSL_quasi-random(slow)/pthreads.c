#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<pthread.h>
#include<gsl/gsl_qrng.h>

struct param {int N,res; gsl_qrng* q;};

void* pi(void* args);

int main(){
	
	int N = 3e7;
	int n1=N/3, n2=N/3, n3=N/3;

	gsl_qrng* qrng1 = gsl_qrng_alloc(gsl_qrng_halton, 2);
	gsl_qrng* qrng2 = gsl_qrng_alloc(gsl_qrng_halton, 2);
	gsl_qrng* qrng3 = gsl_qrng_alloc(gsl_qrng_halton, 2);


	//check to see if the qrngs are different
	//printf("%i\n",*(int*)&qrng1 - *(int*)&qrng2);

	
	struct param p1 = {.N=n1,.res=0,.q=qrng1};
	struct param p2 = {.N=n2,.res=0,.q=qrng2};
	struct param p3 = {.N=n3,.res=0,.q=qrng3};

	pthread_t thread1, thread2, thread3;
	pthread_attr_t* attributes = NULL;

	pthread_create(&thread1, attributes, pi, (void*)&p1);
	pthread_create(&thread2, attributes, pi, (void*)&p2);
	pthread_create(&thread3, attributes, pi, (void*)&p3);


	void* retval = NULL;
	pthread_join(thread1,retval);
	pthread_join(thread2,retval);
	pthread_join(thread3,retval);
		

	gsl_qrng_free(qrng1);
	gsl_qrng_free(qrng2);
	gsl_qrng_free(qrng3);

	//pi((void*)&n1)
	//pi((void*)&n2)
	//pi((void*)&n3)

	double pi_est = (4.0*(p1.res + p2.res + p3.res))/N;
	printf("for N=%i, pi is estimated to be: %.12g\n",N,pi_est);

return 0;
}