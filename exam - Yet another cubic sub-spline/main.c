#include<gsl/gsl_vector.h>
#include"spline.h"
#include<gsl/gsl_interp.h>
#define get gsl_vector_get
#define set gsl_vector_set


void interpolate(double xs[], double ys[], int n, char* path){
	//performs the entire interpolation
	//seperate function so i can do the same analysis for multiple datasets

	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector* p = gsl_vector_alloc(n);
	gsl_vector* c = gsl_vector_alloc(n-1);
	gsl_vector* d = gsl_vector_alloc(n-1);

	FILE* data = fopen(path,"w");


	for (int i=0; i<n; i++){
		set(x,i,xs[i]);
		set(y,i,ys[i]);
		fprintf(data,"%lg %lg\n",xs[i],ys[i]);		//saving datapoints to file
	}
	fprintf(data,"\n\n");

	find_ps(x, y, p);				//finding derivatives at points
	cub_sub(x, y, p, c, d);			//finding cubic sub-spline parameters c and d


	//gsl interp setup for cubic and akima
	gsl_interp* cubic = gsl_interp_alloc(gsl_interp_cspline, n);
	gsl_interp_init(cubic, xs, ys, n);
	gsl_interp* akima = gsl_interp_alloc(gsl_interp_akima, n);
	gsl_interp_init(akima, xs, ys, n);

	
	//perform calculations for points z between given x's
	//used for plotting
	double z = xs[0];
	while (z < xs[n-1]){
		//evaluation, derivative, and integral for sub-spline
		double sz = cub_sub_eval(x, y, p, c, d, z);
		double sz_deriv = cub_sub_deriv(x, y, p, c, d, z);
		double sz_integ = cub_sub_integ(x, y, p, c, d, z);

		//evaluation, derivative, and integral for gsl_cubic
		double cubic_z = gsl_interp_eval(cubic, xs, ys, z,  NULL);
		double cubic_z_deriv = gsl_interp_eval_deriv(cubic, xs, ys, z, NULL);
		double cubic_z_integ = gsl_interp_eval_integ(cubic, xs, ys, xs[0], z, NULL);

		//evaluation, derivative, and integral for akima
		double akima_z = gsl_interp_eval(akima, xs, ys, z,  NULL);
		double akima_z_deriv = gsl_interp_eval_deriv(akima, xs, ys, z, NULL);
		double akima_z_integ = gsl_interp_eval_integ(akima, xs, ys, xs[0], z, NULL);

		//saving all data for point zi to file
		fprintf(data,"%lg %lg %lg %lg ",z, sz, cubic_z, akima_z);
		fprintf(data,"%lg %lg %lg ",sz_deriv, cubic_z_deriv, akima_z_deriv);
		fprintf(data,"%lg %lg %lg\n",sz_integ, cubic_z_integ, akima_z_integ);
		
		z += 0.01;
	}
	fprintf(data,"\n\n");
	fclose(data);

	gsl_interp_free(cubic);
	gsl_interp_free(akima);

	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(p);
	gsl_vector_free(c);
	gsl_vector_free(d);
}


int main(){

	
	int n = 6;
	//example of Akima in action from the book
	//this is to illustrate when regular cubic wiggles too much
	double xs1[] = {-2.5, -1.5, -0.5, 0.5, 1.5, 2.5};
	double ys1[] = {-1, -1, -1, 1, 1, 1};

	interpolate(xs1, ys1, n, "out.data1.txt");

	//example of some random points i chose
	//this is just showing the different interpolations in regular use
	double xs2[] = {0, 1, 2, 3, 4, 5};
	double ys2[] = {0, 2, 5, 8, 5, 7};

	interpolate(xs2, ys2, n, "out.data2.txt");

	
	return 0;
}