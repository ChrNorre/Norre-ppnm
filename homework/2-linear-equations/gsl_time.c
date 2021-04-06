#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_linalg.h>

int main(int argc, char** argv){
	if(argc != 2){
		fprintf(stderr,"not the right amount of arguments to gsl_time (should be 1), is: %d\n",argc-1);
	}
	else {
		int n = atoi(argv[1]);

		gsl_matrix * A = gsl_matrix_alloc(n,n);
		//gsl_matrix * Q = gsl_matrix_alloc(n,n);
		//gsl_matrix * R = gsl_matrix_alloc(n,n);
		gsl_vector * tau = gsl_vector_alloc(n);

		double r = 0;
		for (int i=0; i<n; i++){
			for (int j=0; j<n; j++){
				r = 4.0 * rand() / RAND_MAX;
				gsl_matrix_set(A,j,i,r);
			}
		}
 		gsl_linalg_QR_decomp(A, tau);

		//gsl_linalg_QR_unpack(A, tau, Q, R);

		gsl_matrix_free(A);
		//gsl_matrix_free(Q);
		//gsl_matrix_free(R);
		gsl_vector_free(tau);

		/*
			After QR_decomp, A has technically been factorized, but it is not in a useful form
			I would say that gsl had to provide the Q and R matricies to make a fair comparison, but this would require using QR_unpack
			which requires another two matricies. This makes gsl slower than my qr_decomp, but this is primarialy from the extra alloc
		*/

	}

	return 0;
}