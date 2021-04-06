#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include"qr.h"
#include"tests.h"

void print_matrix(gsl_matrix* C){
	int n = C->size1;
	int m = C->size2;
	for (int i=0; i<n; i++){
		for (int j=0; j<m; j++){
			printf("%10.3g ",gsl_matrix_get(C,i,j));
		}
		printf("\n");
	}
}

void print_vector(gsl_vector* c){
	for (int i=0; i<c->size; i++){
		printf("%10.3g\n",gsl_vector_get(c,i));
	}
}


int main(){
	test_decomp();
	test_solve();
	test_inverse();

return 0;
}