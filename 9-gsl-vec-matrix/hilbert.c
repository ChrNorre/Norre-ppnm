#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>


void vector_print(gsl_vector* v){
	for (int i=0; i < (*v).size; i++){
		printf("%10g ", gsl_vector_get(v,i));
	}
	printf("\n");
}

int main(){
	int n=4;
	gsl_matrix* H = gsl_matrix_alloc(n,n);
	gsl_matrix* Hcopy = gsl_matrix_alloc(n,n);
	gsl_matrix* Hevec = gsl_matrix_alloc(n,n);
	gsl_vector* Heval = gsl_vector_alloc(n);
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);

	//initializing matrix
	for (int i = 0; i < (*H).size1; i++){
		for (int j = 0; j < (*H).size2; j++){
			double Hij = 1.0/(i+j+1); //H is symmetric
			gsl_matrix_set(H,i,j,Hij);
		}
	}
	gsl_matrix_memcpy(Hcopy,H);

	//H is symmetric, so i can use the Real Symmetric functions
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(n);

	gsl_eigen_symmv(Hcopy, Heval, Hevec, w);

	gsl_eigen_symmv_free(w);

	printf("Eigenvectors and Eigenvalues of 4x4 Hilbert matrix\n");
	for (int i = 0; i < (*Heval).size; i++){
		double evali = gsl_vector_get (Heval, i);
		gsl_vector_view eveci = gsl_matrix_column (Hevec, i);
		gsl_vector* x = &eveci.vector;
		printf("\nEigenvector number %i is:\n",i+1);
		vector_print(x);
		printf("With eigenvalue:\t %5g\n",evali);

	
		gsl_blas_dgemv(CblasNoTrans, 1, H, x, 0, y);
		double norm = gsl_blas_dnrm2(y);
		gsl_vector_scale(y, 1.0/norm);
		printf("H*ev%i results in unit vector: \n",i+1);
		vector_print(y);
		printf("With norm of:\t\t %5g\n",norm);		
	}





	
	gsl_matrix_free(H);
	gsl_matrix_free(Hcopy);
	gsl_matrix_free(Hevec);
	gsl_vector_free(Heval);
	gsl_vector_free(x);
	gsl_vector_free(y);
return 0;
}
