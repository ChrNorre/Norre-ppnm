#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include"qr.h"

int main(int argc, char** argv){
	if(argc != 2){
		fprintf(stderr,"not the right amount of arguments to qr_time (should be 1), is: %d\n",argc-1);
	}
	else {
		int n = atoi(argv[1]);

		gsl_matrix * A = gsl_matrix_alloc(n,n);
		gsl_matrix * R = gsl_matrix_alloc(n,n);

		double r = 0;
		for (int i=0; i<n; i++){
			for (int j=0; j<n; j++){
				r = 4.0 * rand() / RAND_MAX;
				gsl_matrix_set(A,j,i,r);
			}
		}
		QR_decomp(A, R);

		gsl_matrix_free(A);
		gsl_matrix_free(R);
	}

	return 0;
}