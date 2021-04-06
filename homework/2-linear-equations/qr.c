#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include"qr.h"


void QR_decomp(gsl_matrix* A, gsl_matrix* R){
	//int n = A->size1;
	int m = A->size2;

	double result;
	for (int i=0; i < m; i++){
		gsl_vector_view ai = gsl_matrix_column(A,i);
		result = gsl_blas_dnrm2(&ai.vector);
		gsl_matrix_set(R,i,i,result);
		gsl_vector_scale(&ai.vector, 1./gsl_matrix_get(R,i,i));
		for (int j=i+1; j < m; j++){
			gsl_vector_view aj = gsl_matrix_column(A,j);
			gsl_blas_ddot(&ai.vector, &aj.vector, &result);
			gsl_matrix_set(R,i,j,result);
			gsl_blas_daxpy(-1*(result), &ai.vector, &aj.vector);
		}
	}
}

void QR_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	//int n = Q->size1;
	int m = Q->size2;

	gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);   //x = 1*QT*b + 0*x

	//back-propagation
	for(int i=m-1; i>=0; i--){
		double s = gsl_vector_get(x,i);
		for (int k=i+1; k<m; k++){
			s -= gsl_matrix_get(R, i, k) * gsl_vector_get(x, k);			
		}
		gsl_vector_set(x, i, s/gsl_matrix_get(R, i, i)); 
	}

}

void QR_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B, gsl_vector* x){	
	int n = Q->size1;

	for (int i = 0; i < n; i++){
		gsl_vector_view rowi = gsl_matrix_column(B,i);
		QR_solve(Q, R, &rowi.vector, x);
		gsl_vector_memcpy(&rowi.vector, x);
	}
}