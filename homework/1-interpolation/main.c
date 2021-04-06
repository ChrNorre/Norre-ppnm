#include<stdio.h>
#include<math.h>
#include<gsl/gsl_interp.h>
#include<assert.h>
#include<gsl/gsl_vector.h>


int binsearch(int n, double* x, double z){/* locates the interval for z by bisection */ 
	assert(x[0]<=z && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
	}
	return i;
}

int binsearch_vec(gsl_vector* x, double z){/* locates the interval for z by bisection */ 
	assert(gsl_vector_get(x,0)<=z && z<=gsl_vector_get(x,x->size-1));
	int i=0, j=x->size-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>gsl_vector_get(x,mid)) i=mid; else j=mid;
	}
	return i;
}

//implemented in linear.c
double linterp(int n, double* x, double* y, double z);
double linterp_integ(int n, double* x, double* y, double z);


//implemented in quadratic.c
void quaterp(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c);
double quaterp_eval(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c, double z);
double quaterp_integ(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c, double z);
double quaterp_diff(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c, double z);


//implemented in cubic.c
void cubterp(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c,gsl_vector *d);
double cubterp_eval(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c,gsl_vector *d, double z);
double cubterp_integ(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c,gsl_vector *d, double z);
double cubterp_diff(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c,gsl_vector *d, double z);


int main(){

	int n = 11;
	double a=-2, bint=2;
	double xa[n], ya[n];
	gsl_vector *x = gsl_vector_alloc(n);
	gsl_vector *y = gsl_vector_alloc(n);

	

	//making x and y data
	FILE* outfile = fopen("datapoints.txt","w");
	for(int i=0; i<n; i++){
		xa[i] = a+(bint-a)*i/(n-1);
		ya[i] = 1/(1+pow(xa[i],2));
		gsl_vector_set(x,i,xa[i]);
		gsl_vector_set(y,i,ya[i]);
		fprintf(outfile, "%g %g\n",xa[i],ya[i]);
	}
	fclose(outfile);	



	//linear interp
	gsl_interp* linear = gsl_interp_alloc(gsl_interp_linear, n);
	gsl_interp_init(linear, xa, ya, n);

	double z, myz, gslz, myzint, gslint;
	FILE* outlinear = fopen("linear.txt","w");
	for (int i=0; i < 1000; i++){
		z = -2 + 4*(1.0*i/1000);
		myz = linterp(n,xa,ya,z);
		gslz = gsl_interp_eval(linear, xa, ya, z,  NULL);
		myzint = linterp_integ(n, xa, ya, z);
		gslint = gsl_interp_eval_integ(linear, xa, ya, xa[0], z, NULL);
		fprintf(outlinear,"%g %g %g %g %g\n",z,myz,gslz,myzint,gslint);
	}
	fclose(outlinear);
	gsl_interp_free(linear);





	//quadratic interp
	double myzy_qua, myzint_qua, myzdiff_qua;
	gsl_vector* b_qua = gsl_vector_alloc(n-1);
	gsl_vector* c_qua = gsl_vector_alloc(n-1);
	FILE* outquadratic = fopen("quadratic.txt","w");
	for (int i=0; i < 1000; i++){
		z = -2 + 4*(1.0*i/1000);
		quaterp(x,y,b_qua,c_qua);
		myzy_qua = quaterp_eval(x,y,b_qua,c_qua,z);
		myzint_qua = quaterp_integ(x,y,b_qua,c_qua,z);
		myzdiff_qua = quaterp_diff(x,y,b_qua,c_qua,z);
		fprintf(outquadratic,"%g %g %g %g\n",z,myzy_qua,myzint_qua,myzdiff_qua);
	}
	fclose(outquadratic);
	gsl_vector_free(b_qua);
	gsl_vector_free(c_qua);






	//cubic interp
	gsl_interp* cubic = gsl_interp_alloc(gsl_interp_cspline, n);
	gsl_interp_init(cubic, xa, ya, n);

	gsl_vector* b = gsl_vector_alloc(n);
	gsl_vector* c = gsl_vector_alloc(n-1);
	gsl_vector* d = gsl_vector_alloc(n-1);

	double myz_cub, gslz_cub, myzint_cub, gslint_cub, myzdiff_cub, gsldiff_cub;
	FILE* outcubic = fopen("cubic.txt","w");
	for (int i=0; i < 1000; i++){
		z = -2 + 4*(1.0*i/1000);
		cubterp(x,y,b,c,d);
		myz_cub = cubterp_eval(x,y,b,c,d,z);
		gslz_cub = gsl_interp_eval(cubic, xa, ya, z,  NULL);
		myzint_cub = cubterp_integ(x,y,b,c,d,z);
		gslint_cub = gsl_interp_eval_integ(cubic, xa, ya, xa[0], z, NULL);
		myzdiff_cub = cubterp_diff(x,y,b,c,d,z);
		gsldiff_cub = gsl_interp_eval_deriv(cubic, xa, ya, z, NULL);
		fprintf(outcubic,"%g %g %g %g %g %g %g\n",z,myz_cub,gslz_cub,myzint_cub,gslint_cub,myzdiff_cub,gsldiff_cub);
	}
	fclose(outcubic);
	gsl_interp_free(cubic);
	gsl_vector_free(b);
	gsl_vector_free(c);
	gsl_vector_free(d);




	gsl_vector_free(x);
	gsl_vector_free(y);

	return 0;
}