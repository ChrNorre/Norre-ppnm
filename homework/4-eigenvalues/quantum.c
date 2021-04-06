#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include"jacobi.h"


int main(){
	//number of points on line
	int n=20;

	//initialize hamiltonian
	double s=1.0/(n+1);
	gsl_matrix* H = gsl_matrix_alloc(n,n);
	for(int i=0;i<n-1;i++){
		gsl_matrix_set(H,i,i,-2);
		gsl_matrix_set(H,i,i+1,1);
		gsl_matrix_set(H,i+1,i,1);
	}
	gsl_matrix_set(H,n-1,n-1,-2);
	gsl_matrix_scale(H,-1/s/s);

	//diagonaize hamiltonian
	gsl_matrix* V = gsl_matrix_alloc(n,n);
	jacobi_diag(H,V);


	//checking energies
	printf("\nEigenenergies: num, caluclated, exact\n");
	for (int k=0; k < n/3; k++){
		double exact = M_PI*M_PI*(k+1)*(k+1);
		double calculated = gsl_matrix_get(H,k,k);
		printf("%i %7g %7g\n",k,calculated,exact);
	}
	printf("\n\n");

	FILE* outdata = fopen("data.txt","w");	
	for (int i = 0; i < n; i++){
		fprintf(outdata,"0 0 0 0 0 0\n");
		fprintf(outdata,"%10g",(i+1.0)/(n+1));
		for (int k = 0; k < 4; k++){
			fprintf(outdata," %10g",gsl_matrix_get(V,i,k));
		}
		fprintf(outdata,"\n1 0 0 0 0 0\n");
	}
		

	gsl_matrix_free(H);
	gsl_matrix_free(V);


	return 0;
}