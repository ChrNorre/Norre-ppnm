#include<math.h>
#include<gsl/gsl_qrng.h>

struct param {int N,res; gsl_qrng* q;};

void* pi(void* args){
	struct param * p= (struct param *) args;

	gsl_qrng * qrng = (*p).q;

	double x[2];
	int ncirc = 0;
	for(int i = 0; i<(*p).N; i++){
		gsl_qrng_get(qrng, x);
		if (pow(x[0],2) + pow(x[1],2) <= 1){
			ncirc++;
		}
	}
	(*p).res = ncirc;

return NULL;
}