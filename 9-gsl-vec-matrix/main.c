#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>


void vector_print(char* s, gsl_vector* v){
	printf("%s\n",s);
	for (int i=0; i < (*v).size; i++){
		printf("%10g ", gsl_vector_get(v,i));
	}
	printf("\n");
}

int main(){
	int n=3;
	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* Acopy = gsl_matrix_alloc(n,n);
	gsl_vector* b = gsl_vector_alloc(n);
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);

	//initializing matrix
	double lsA[9] = {6.13, -2.90, 5.86, 8.08, -6.31, -3.89, -4.36, 1.00, 0.19};
	for (int i = 0; i < (*A).size1; i++){
		for (int j = 0; j < (*A).size2; j++){
			gsl_matrix_set(A,i,j,lsA[j+3*i]);
		}
	}
	gsl_matrix_memcpy(Acopy,A);

	//initializing vector
	double lsb[3] = {6.23, 5.37, 2.29};
	for (int i = 0; i < (*b).size; i++){
		gsl_vector_set(b,i,lsb[i]);
	}


	//solving Ax=b
	gsl_linalg_HH_solve(Acopy,b,x);

	//y = 1*A*x + 0*y          here y = A*x	
	gsl_blas_dgemv(CblasNoTrans, 1, A, x, 0, y);


	//if x is correct, then these should be the same
	vector_print("b=",b);
	vector_print("Ax=",y);



	gsl_matrix_free(A);
	gsl_matrix_free(Acopy);
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_vector_free(y);
return 0;
}
