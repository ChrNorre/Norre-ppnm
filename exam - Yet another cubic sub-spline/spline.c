#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include"spline.h"
#define get gsl_vector_get
#define set gsl_vector_set


int binsearch_vec(gsl_vector* x, double z){
	// locates the interval for z by bisection 
	int i=0, j=x->size-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>gsl_vector_get(x,mid)) i=mid; else j=mid;
	}
	return i;
}

void find_ps(gsl_vector* x, gsl_vector* y, gsl_vector* p){
	// x and y are the given data, p is empty vector which will get filled with the derivatives
	// using a quadratic spline
	// i am not using this spline, i was only interested in the derivatives p

	int n = x->size;

	double xm, ym, xp, yp, denum, bi, ci;
	xp = get(x,1)-get(x,0);
	yp = get(y,1)-get(y,0);

	for (int i = 1; i < n-1; i++){
		xm = -xp;
		ym = -yp;
		xp = get(x,i+1)-get(x,i);
		yp = get(y,i+1)-get(y,i);

		denum = xm*xp*(xm-xp);
		bi =  (yp*xm*xm - ym*xp*xp) / denum;
		ci = -(yp*xm    - ym*xp   ) / denum;		//only matters for c1 and c(n-2)
		
		if (i==1){
			set(p,0,bi + 2*ci*xm);
		}
		if (i==n-2){
			set(p,n-1,bi + 2*ci*xp);
		}
		set(p,i,bi);
	}
}

void cub_sub(gsl_vector* x, gsl_vector* y, gsl_vector* p, gsl_vector* c, gsl_vector* d){
	double xp, yp, ci, di;
	for (int i=0; i < c->size; i++){
		xp = get(x,i+1) - get(x,i);
		yp = get(y,i+1) - get(y,i);
		ci = (3*yp - (get(p,i+1) + 2*get(p,i))*xp) / (xp*xp);
		di = ((get(p,i+1) + get(p,i))*xp - 2*yp) / (xp*xp);
		set(c,i,ci);
		set(d,i,di);
	}
}

double cub_sub_eval(gsl_vector* x, gsl_vector* y, gsl_vector* p, gsl_vector* c, gsl_vector* d, double z){
	int i = binsearch_vec(x,z);		//binary search to find which interval z lies in
	double zxi = z - gsl_vector_get(x,i);
	double Siz = gsl_vector_get(y,i) + zxi*(gsl_vector_get(p,i) + zxi*(gsl_vector_get(c,i) + zxi*gsl_vector_get(d,i)));
	return Siz;
}


double cub_sub_deriv(gsl_vector* x, gsl_vector* y, gsl_vector* p, gsl_vector* c, gsl_vector* d, double z){
	int i = binsearch_vec(x,z);
	double xp = z - get(x,i);
	return get(p,i) + 2*get(c,i)*xp + 3*get(d,i)*xp*xp;
}


double cub_sub_integ2(gsl_vector* x, gsl_vector* y, gsl_vector* p, gsl_vector* c, gsl_vector* d, double z){
	int i = binsearch_vec(x,z);
	double res = 0;
	double xp;
	for (int j=0; j<i; j++){
		xp = get(x,j+1) - get(x,j);
		res += 1.0/12 * xp * ( 3*get(y,j+1) + 9*get(y,j) + xp*(3*get(p,j) + xp*get(c,j)));
	}
	double S_z = cub_sub_eval(x, y, p, c, d, z);
	xp = z - get(x,i);
	res += 1.0/12 * xp * ( 3*S_z + 9*get(y,i) + xp*(3*get(p,i) + xp*get(c,i)));
	return res;
}

double cub_sub_integ(gsl_vector* x, gsl_vector* y, gsl_vector* p, gsl_vector* c, gsl_vector* d, double z){
	int i = binsearch_vec(x,z);
	double res = 0;
	double xp;
	for (int j=0; j<i; j++){
		xp = get(x,j+1) - get(x,j);
		res += 1.0/12 * xp * ( 12*get(y,j) + xp*(6*get(p,j) + xp*(4*get(c,j) + 3*xp*get(d,j))));
	}
	xp = z - get(x,i);
	res += 1.0/12 * xp * ( 12*get(y,i) + xp*(6*get(p,i) + xp*(4*get(c,i) + 3*xp*get(d,i))));
	return res;
}