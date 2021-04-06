#include<gsl/gsl_vector.h>

int binsearch_vec(gsl_vector* x, double z);

void cubterp(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c,gsl_vector *d){
	int n = x->size;
	int i;
	double  p[n-1], h[n-1];
	for (i=0; i < n-1; i++){
		h[i] =  gsl_vector_get(x,i+1) - gsl_vector_get(x,i);
		p[i] = (gsl_vector_get(y,i+1) - gsl_vector_get(y,i)) / h[i];
	}

	double D[n] , Q[n-1] , B[n];
	D[0] = 2;
	Q[0] = 1;
	B[0] = 3*p[0];
	for (i=0; i < n-2; i++){
		D[i+1] = 2 * h[i]/h[i+1] + 2;
		Q[i+1] = h[i]/h[i+1];
		B[i+1] = 3 * (p[i] + p[i+1] * h[i]/h[i+1]);
	}
	D[x->size-1] = 2;
	B[x->size-1] = 3*p[n-2];


	for (i = 1; i < n; i++){
		D[i] -= Q[i-1] / D[i-1];
		B[i] -= B[i-1] / D[i-1];
	}

	double bi;
	gsl_vector_set(b,n-1, B[n-1] / D[n-1]);
	for (i = n-2; i >= 0; i--) {
		bi = (B[i] - Q[i] * gsl_vector_get(b,i+1)) / D[i];
		gsl_vector_set(b, i, bi);
	}

	double bi1;
	for (i = 0; i < n-1; i++){
		bi = gsl_vector_get(b,i);
		bi1 = gsl_vector_get(b,i+1);
		gsl_vector_set(c,i, (-2*bi - bi1 + 3*p[i])/h[i] );
		gsl_vector_set(d,i, (bi + bi1 - 2*p[i])/(h[i]*h[i]) );
	}
}

double cubterp_eval(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c,gsl_vector *d, double z){
	int i = binsearch_vec(x,z);
	double zxi = z - gsl_vector_get(x,i);
	double siz = gsl_vector_get(y,i) + zxi*(gsl_vector_get(b,i) + zxi*(gsl_vector_get(c,i) + gsl_vector_get(d,i)*zxi));
	return siz;
}


double cubterp_integ(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c,gsl_vector *d, double z){
	int i = binsearch_vec(x,z);
	double res = 0;
	double delx;
	for (int j=0; j < i; j++){
		delx = gsl_vector_get(x,j+1) - gsl_vector_get(x,j);
		res += 1./12 * delx*( 3*gsl_vector_get(y,j+1) + 9*gsl_vector_get(y,j) + 3*delx*gsl_vector_get(b,j) + delx*delx*gsl_vector_get(c,j) );
	}
	double zy = cubterp_eval(x,y,b,c,d,z);
	delx = z - gsl_vector_get(x,i);
	res += 1./12 * delx*( 3*zy + 9*gsl_vector_get(y,i) + 3*delx*gsl_vector_get(b,i) + delx*delx*gsl_vector_get(c,i) );
	return res;
}

double cubterp_diff(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c,gsl_vector *d, double z){
	int i = binsearch_vec(x,z);
	double delx = z - gsl_vector_get(x,i);
	double res = gsl_vector_get(b,i) + 2*gsl_vector_get(c,i)*delx + 3*gsl_vector_get(d,i)*delx*delx;
	return res;
}