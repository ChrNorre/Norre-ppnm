#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include"integ.h"
#include<stdio.h>


double integrate(double (*f)(double, int*), double a, double b, double abs, double rel, int* counter){
	/* Recursive adaptive integrator for open quadrature
	Abscissas: for [0,1], {1/6, 2/6, 4/6, 5/6}
	4-weights: {2/6, 1/6, 1/6, 2/6}
	2-weights: {1/4, 1/4, 1/4, 1/4}
	*/

	double p26 = f(a + (b-a)*(2.0/6), counter);
	double p46 = f(a + (b-a)*(4.0/6), counter);

	return integrate24(f, a, b, abs, rel, p26, p46, counter);
}

double integrate24(double (*f)(double, int*), double a, double b, double abs, double rel, double p26, double p46, int* counter){

	double p16 = f(a + (b-a)*(1.0/6), counter);
	double p56 = f(a + (b-a)*(5.0/6), counter);

	double Q = (b-a)*(2*p16 + 1*p26 + 1*p46 + 2*p56)/6; //sum_i   w_i * f(x_i)
	double q = (b-a)*(1*p16 + 1*p26 + 1*p46 + 1*p56)/4; //sum_i   w_alternative_i * f(x_i)

	double tol = abs + rel*fabs(Q);
	double error = fabs(Q-q);

	if (error < tol){
		return Q;
	} 
	else {
		double Qleft = integrate24(f, a, (a+b)/2, abs/sqrt(2.0), rel, p16, p26, counter);
		double Qright = integrate24(f, (a+b)/2, b, abs/sqrt(2.0), rel, p46, p56, counter);
		return Qleft+Qright;
	}

}







//this is a very inelegant solution but it is the only thing that works when i want to
// pass a function as a parameter
double a_field;
double b_field;
double (*f_field)(double, int*);

//variable transform function f into Clenshaw-Curtis	
double CC(double x, int* counter){
	double newX = (a_field+b_field)/2 + (b_field-a_field)/2 * cos(x);
	double outside = sin(x)*(b_field-a_field)/2;
	return f_field(newX, counter)*outside;
}

double integrate_CC(double (*f)(double, int*), double a, double b, double abs, double rel, int* counter){
	/* Recursive adaptive integrator for open quadrature
	Clenshaw-Curtis variable transformation
	*/
	
	a_field = a;
	b_field = b;
	f_field = f;

	
	double p26 = CC(0 + (M_PI-0)*(2.0/6), counter);
	double p46 = CC(0 + (M_PI-0)*(4.0/6), counter);

	//pass CC and new boundary [a,b] -> [0,pi]
	return integrate24_CC(&CC, 0, M_PI, abs, rel, p26, p46, counter);
}

double integrate24_CC(double (*CC)(double, int*), double a, double b, double abs, double rel, double p26, double p46, int* counter){

	
	double p16 = CC(a + (b-a)*(1.0/6), counter);
	double p56 = CC(a + (b-a)*(5.0/6), counter);
	
	double Q = (b-a)*(2*p16 + 1*p26 + 1*p46 + 2*p56)/6; //sum_i   w_i * f(x_i)
	double q = (b-a)*(1*p16 + 1*p26 + 1*p46 + 1*p56)/4; //sum_i   w_alternative_i * f(x_i)

	double tol = abs + rel*fabs(Q);
	double error = fabs(Q-q);

	if ((error < tol) || ((b-a) < 1e-6)) {
		// I added the (b-a) < 1e-6 condition because i was getting segmentation faults without it
		// or -inf if (b-a) < ~1e-8 or smaller
		//this was strangely only a problem for log(x)/sqrt(x) and not for 1/sqrt(x)
		return Q;
	} 
	else {
		double Qleft = integrate24_CC(CC, a, (a+b)/2, abs/sqrt(2.0), rel, p16, p26, counter);
		double Qright = integrate24_CC(CC, (a+b)/2, b, abs/sqrt(2.0), rel, p46, p56, counter);
		return Qleft+Qright;
	}

}