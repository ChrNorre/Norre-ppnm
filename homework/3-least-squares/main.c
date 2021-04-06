#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include"least_squares_fit.h"


double funs(int i, double x){
	switch(i){
		case 0: return 1; break;
		case 1: return x; break;
		case 2: return x*x; break;
		case 3: return x*x*x; break;
		default: return NAN;
	}
}

int main(){

	double t[9] = {1, 2, 3, 4, 6, 9, 10, 13, 15};
	double activity[9] = {117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1};


	int n = 9;
	int k = 2;


	gsl_vector* xs = gsl_vector_alloc(n);
	gsl_vector* ys = gsl_vector_alloc(n);
	gsl_vector* dys = gsl_vector_alloc(n);
	gsl_vector* cs = gsl_vector_alloc(k);
	gsl_matrix* Sigma = gsl_matrix_alloc(k,k);


	for (int i = 0; i < n; i++){
		gsl_vector_set(xs,i,t[i]);
		gsl_vector_set(ys,i,log(activity[i]));
		gsl_vector_set(dys,i,1.0/20);
	}


	least_squares(xs, ys, dys, &funs, cs, Sigma);

	printf("lambda = %g +/- %g\n", -1*gsl_vector_get(cs,1), sqrt(gsl_matrix_get(Sigma,1,1)));
	printf("half-life = ln(2)/lambda = from %g to %g (days)\n", -1*log(2)/(gsl_vector_get(cs,1)-sqrt(gsl_matrix_get(Sigma,1,1))), -1*log(2)/(gsl_vector_get(cs,1)+sqrt(gsl_matrix_get(Sigma,1,1))));
	printf("half-life of radium-224 from wikipedia = 3.63 (days)\n");


	FILE* outdata = fopen("data.txt","w");
	for (int i=0; i < n; i++){
		fprintf(outdata,"%10g %10g %10g\n",gsl_vector_get(xs,i), gsl_vector_get(ys,i), gsl_vector_get(dys,i));
	}
	fclose(outdata);


	int m = 200;
	double x, y, ym, yp;

	gsl_vector* csm = gsl_vector_alloc(k);
	gsl_vector* csp = gsl_vector_alloc(k);
	gsl_vector_set(csm, 0, gsl_vector_get(cs,0));
	gsl_vector_set(csp, 0, gsl_vector_get(cs,0));
	gsl_vector_set(csm, 1, gsl_vector_get(cs,1) - gsl_matrix_get(Sigma, 1, 1));
	gsl_vector_set(csp, 1, gsl_vector_get(cs,1) + gsl_matrix_get(Sigma, 1, 1));

	FILE* outfit = fopen("fit.txt","w");
	for (int i=0; i < m; i++){
		x = 1 + 14* (1.0*i/m);
		y = eval_fun(&funs, cs, x);
		ym = eval_fun(&funs, csm, x);
		yp = eval_fun(&funs, csp, x);
		fprintf(outfit,"%10g %10g %10g %10g\n", x, y, ym, yp);
	}
	fclose(outfit);

	gsl_vector_free(csm);
	gsl_vector_free(csp);


	gsl_vector_free(xs);
	gsl_vector_free(ys);
	gsl_vector_free(dys);
	gsl_vector_free(cs);
	gsl_matrix_free(Sigma);


	return 0;
}