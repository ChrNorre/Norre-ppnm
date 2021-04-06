#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include"qr.h"


void print_vector(gsl_vector* c);

void print_matrix(gsl_matrix* C);

void test_decomp(){

	int n = 5;
	int m = 3;

	gsl_matrix * A = gsl_matrix_alloc(n,m);
	gsl_matrix * Q = gsl_matrix_alloc(n,m);
	gsl_matrix * R = gsl_matrix_alloc(m,m);

	double r = 0;
	for (int i=0; i<m; i++){
		for (int j=0; j<n; j++){
			r = 4.0 * rand() / RAND_MAX;
			gsl_matrix_set(A,j,i,r);
			gsl_matrix_set(Q,j,i,r);
		}
	}
	QR_decomp(Q, R);	//Q=A

	printf("\nTesting for QR-decomposition:\n");
	printf("R matrix = \n");
	print_matrix(R);

	printf("\nA = \n");
	print_matrix(A);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Q, R, 0, A);
	printf("\nQ*R = \n");
	print_matrix(A);

	//gsl_matrix_memcpy(A, Q);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, Q, Q, 0, R);
	printf("\nQT * Q = \n");
	print_matrix(R);
	printf("\n\n");


	gsl_matrix_free(A);
	gsl_matrix_free(Q);
	gsl_matrix_free(R);

}


void test_solve(){

	int n = 4;
	int m = 4;

	gsl_matrix * A = gsl_matrix_alloc(n,m);
	gsl_matrix * Q = gsl_matrix_alloc(n,m);
	gsl_matrix * R = gsl_matrix_alloc(m,m);
	gsl_vector * x = gsl_vector_alloc(m);
	gsl_vector * b = gsl_vector_alloc(n);

	double r = 0;
	for (int i=0; i<n; i++){
		for (int j=0; j<m; j++){
			r = 4.0 * rand() / RAND_MAX;
			gsl_matrix_set(A,i,j,r);
			gsl_matrix_set(Q,i,j,r);
		}
		r = 4.0 * rand() / RAND_MAX;
		gsl_vector_set(b,i,r);
	}

	QR_decomp(Q, R);

	printf("\nTesting for QR-Solve:\n");
	printf("A matrix = \n");
	print_matrix(A);

	printf("\nb vector = \n");
	print_vector(b);


	QR_solve(Q,R,b,x);

	printf("\nQRx=b    =>     x = \n");
	print_vector(x);


	gsl_blas_dgemv(CblasNoTrans, 1, A, x, 0, b);
	printf("\nA*x =      (should be b)\n");
	print_vector(b);
	printf("\n\n");
	

	gsl_matrix_free(A);
	gsl_matrix_free(Q);
	gsl_matrix_free(R);

	gsl_vector_free(x);
	gsl_vector_free(b);

}

void test_inverse(){

	int n = 5;

	gsl_matrix * A = gsl_matrix_alloc(n,n);
	gsl_matrix * Q = gsl_matrix_alloc(n,n);
	gsl_matrix * R = gsl_matrix_alloc(n,n);
	gsl_matrix * B = gsl_matrix_calloc(n,n);

	gsl_vector* x = gsl_vector_alloc(n);

	double r = 0;
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			r = 4.0 * rand() / RAND_MAX;
			gsl_matrix_set(A,i,j,r);
			gsl_matrix_set(Q,i,j,r);
		}
		r = 4.0 * rand() / RAND_MAX;
		gsl_matrix_set(B,i,i,1);
	}

	QR_decomp(Q, R);
	QR_inverse(Q, R, B, x);

	printf("\nTesting for QR-Inverse:\n");
	printf("A matrix = \n");
	print_matrix(A);

	printf("\nA-1 matrix = \n");
	print_matrix(B);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, B, 0, Q);

	printf("\nA*A-1 =      (should be I)\n");
	print_matrix(Q);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, B, A, 0, Q);

	printf("\nA-1*A =      (should be I)\n");
	print_matrix(Q);
	printf("\n\n");
	

	gsl_matrix_free(A);
	gsl_matrix_free(Q);
	gsl_matrix_free(R);

	gsl_vector_free(x);

}