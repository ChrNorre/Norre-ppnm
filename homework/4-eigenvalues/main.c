#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include"jacobi.h"

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
	int n = 4;

	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* V = gsl_matrix_alloc(n,n);
	gsl_matrix* A_copy = gsl_matrix_alloc(n,n);
	gsl_matrix* C = gsl_matrix_alloc(n,n);
	gsl_matrix* B = gsl_matrix_alloc(n,n);
	

	//V unit matrix, A random real symmetric
	//double r = 1.0*rand()/RAND_MAX;
	double r;
	for (int i = 0; i < n; i++){
		for (int j = i; j < n; j++){
			r = 1.0*rand()/RAND_MAX;
			if (i==j){
				gsl_matrix_set(A,i,j,r);
			}
			else {
				gsl_matrix_set(A,i,j,r);
				gsl_matrix_set(A,j,i,r);
			}
		}
	}

	gsl_matrix_memcpy(A_copy, A);

	printf("\nMatrix A:\n");
	print_matrix(A);

	jacobi_diag(A, V);


	//checking V * D * V^T = A
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, A, V, 0, C);	//C = D * V^T
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, V, C, 0, B); // B = V * C = V * D * V^T
	printf("\nV * D * V^T (should be A):\n");
	print_matrix(B);


	printf("\nMatrix D (should be diagonal):\n");
	print_matrix(A);

	//checking V^T * A * V = D
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A_copy, V, 0, C);	//C = A * V
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, C, 0, B); // B = V^T * C = V^T * A * V
	printf("\nV^T * A * V (should be D):\n");
	print_matrix(B);


	//checking V^T * V = 1
	gsl_matrix_memcpy(B, V);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, B, V, 0, C);	//C = V^T * V
	printf("\nV^T * V (should be unit):\n");
	print_matrix(C);
	printf("\n\n");



	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_matrix_free(A_copy);
	gsl_matrix_free(C);
	gsl_matrix_free(B);

	return 0;
}