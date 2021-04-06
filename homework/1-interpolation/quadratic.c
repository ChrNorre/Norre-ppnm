#include<gsl/gsl_vector.h>
#include<stdio.h>

int binsearch_vec(gsl_vector* x, double z);

void quaterp(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c){
	int i;
	double ci;
	double  p[c->size], h[c->size];
	for (i=0; i < x->size-1; i++){
		h[i] = gsl_vector_get(x,i+1) - gsl_vector_get(x,i);
		p[i] = (gsl_vector_get(y,i+1) - gsl_vector_get(y,i)) / h[i];
	}

	gsl_vector_set(c,0,0);
	for (i=0; i < c->size-1; i++){
		ci = (p[i+1] - p[i] - gsl_vector_get(c,i) * h[i]) / h[i+1];
		gsl_vector_set(c, i+1, ci);
	}

	gsl_vector_set(c,c->size-1, gsl_vector_get(c, c->size-1) / 2);
	for (i=c->size-2; i >= 0; i--){
		ci = (p[i+1] - p[i] - gsl_vector_get(c,i+1) * h[i+1]) / h[i];
		gsl_vector_set(c, i, ci);
	}

	for (i = 0; i < c->size; i++){
		gsl_vector_set(b,i, p[i] - gsl_vector_get(c,i)*h[i]);
	}
}

double quaterp_eval(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c, double z){
	int i = binsearch_vec(x,z);
	double zxi = z - gsl_vector_get(x,i);
	double siz = gsl_vector_get(y,i) + gsl_vector_get(b,i)*zxi + gsl_vector_get(c,i)*zxi*zxi;
	return siz;
}


double quaterp_integ(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c, double z){
	int i = binsearch_vec(x,z);
	double res = 0;
	double delx;
	for (int j=0; j < i; j++){
		delx = gsl_vector_get(x,j+1) - gsl_vector_get(x,j);
		res += 1./6 * delx * ( 2 * gsl_vector_get(y,j+1) + 4 * gsl_vector_get(y,j) + delx * gsl_vector_get(b,j) );
	}
	double zy = quaterp_eval(x,y,b,c,z);
	delx = z - gsl_vector_get(x,i);
	res += 1./6 * delx * ( 2 * zy + 4 * gsl_vector_get(y,i) + delx * gsl_vector_get(b,i) );
	return res;
}

double quaterp_diff(gsl_vector *x,gsl_vector *y,gsl_vector *b,gsl_vector *c, double z){
	int i = binsearch_vec(x,z);
	double delx = z - gsl_vector_get(x,i);
	double res = gsl_vector_get(b,i) + 2*gsl_vector_get(c,i) * delx;
	return res;
}