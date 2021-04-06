#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include"jacobi.h"

int equal(double a, double b, double tau, double epsilon){
	if (fabs(a-b) < tau) return 1;
	else if (fabs(a-b)/(fabs(a)+fabs(b)) < epsilon/2) return 1;
	else return 0;
}

void timesJ(gsl_matrix* A, int p, int q, double theta){
	// A = A * J(p,q,theta)

	double s = sin(theta);
	double c = cos(theta);
	for (int i = 0; i < A->size1; i++){
		double aip = c*gsl_matrix_get(A,i,p) - s*gsl_matrix_get(A,i,q);
		double aiq = s*gsl_matrix_get(A,i,p) + c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,aip);
		gsl_matrix_set(A,i,q,aiq);
	}
}

void Jtimes(gsl_matrix* A, int p, int q, double theta){
	// A = J(p,q,theta) * A

	double s = sin(theta);
	double c = cos(theta);
	for (int i = 0; i < A->size2; i++){
		double api = c*gsl_matrix_get(A,p,i) + s*gsl_matrix_get(A,q,i);
		double aqi = -s*gsl_matrix_get(A,p,i) + c*gsl_matrix_get(A,q,i);
		gsl_matrix_set(A,p,i,api);
		gsl_matrix_set(A,q,i,aqi);
	}

}

void jacobi_diag(gsl_matrix* A, gsl_matrix* V){
	// makes several sweeps to diagonalize A (real symmetric)
	// Eigenvectors will be stored in V

	//initializing V to unit matrix
	for (int i = 0; i < V->size1; i++){
		for (int j = i; j < V->size2; j++){
			if (i==j){
				gsl_matrix_set(V,i,j,1);
			}
			else {
				gsl_matrix_set(V,i,j,0);
				gsl_matrix_set(V,j,i,0);
			}
		}
	}
	
	//diagonal elements become eigenvalues, stop when a sweep does not change the diagonal values
	int changed;
	do {
		changed = 0;
		//single sweep, run over upper triangular part
		for (int p = 0; p < A->size1-1; p++){
			for (int q = p+1; q < A->size1; q++){
				//find theta
				double apq = gsl_matrix_get(A,p,q);
				double app = gsl_matrix_get(A,p,p);
				double aqq = gsl_matrix_get(A,q,q);
				double theta = atan2(2*apq, aqq - app)/2;

				//manual calculation of new diagonal elements to see if rotation is necessary
				//alternatively, i could rotate and see if anything changed, but that is wasted O(n^3) work for the last sweep
				//compared to an extra O(1) work per sweep to check if it is even needed.
				double c = cos(theta);
				double s = sin(theta);
				double new_app = c*c*app - 2*s*c*apq + s*s*aqq;
				double new_aqq = s*s*app + 2*s*c*apq + c*c*aqq;

				//if(equal(new_app, app, 1.0e-16, 1.0e-16)==0 || equal(new_aqq, aqq, 1.e-16,  1.e-16)==0) {
				if (new_app != app || new_aqq != aqq){
					//printf("yes\n");
					changed = 1;
					// do A = J^T * A * J (for this p and q)
					timesJ(A,p,q,theta);
					Jtimes(A,p,q,-theta); //minus theta because transposed

					//calculate V (this is the entire reason for seperating timesJ and Jtimes)
					// V = V * J
					timesJ(V,p,q,theta);
				}
			}
		}
	} while(changed != 0);
}