#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<omp.h>
#include<pthread.h>	//this is needed to make time(NULL) work

struct param {int N,res,seednum;};

void* pi(void* args){
	struct param * p= (struct param *) args;

	unsigned int seed1 = (unsigned) time(NULL) + (*p).seednum + 1;
	unsigned int seed2 = (unsigned) time(NULL) + (*p).seednum + 2;
	//printf("%i\n%i\n",seed1,seed2);

	double x1,x2;
	int ncirc = 0;
	for(int i = 0; i<(*p).N; i++){
		x1 = (1.0*rand_r(&seed1))/RAND_MAX;
		x2 = (1.0*rand_r(&seed2))/RAND_MAX;
		if (pow(x1,2) + pow(x2,2) <= 1){
			ncirc++;
		}
	}
	(*p).res = ncirc;

return NULL;
}


int main(){
	
	int N = 3e8;
	int n1=N/3, n2=N/3, n3=N/3;

	struct param p1 = {.N=n1,.res=0,.seednum=10};
	struct param p2 = {.N=n2,.res=0,.seednum=20};
	struct param p3 = {.N=n3,.res=0,.seednum=30};


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
	

	double pi_est = (4.0*(p1.res + p2.res + p3.res))/N;
	printf("for N=%i, pi is estimated to be: %.12g\n",N,pi_est);

return 0;
}