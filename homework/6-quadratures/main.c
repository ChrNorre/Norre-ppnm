#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include"integ.h"
#include<stdio.h>
#include<gsl/gsl_integration.h>




int main(){


	double a = 0.0;
	double b = 1.0;

	
	double f_test(double x, int* counter){
		(*counter)++;
		return sqrt(x);
	}
	double f_test2(double x, int* counter){
		(*counter)++;
		return 4*sqrt(1-x*x);
	}

	int counter = 0;
	double err = 0;	
	double abs = 1e-10;
	double rel = 1e-10;

	double res = integrate(&f_test, a, b, abs, rel, &counter, &err);
	printf("\nint(sqrt(x))[0,1] = 2/3 = \t %lg +/- %lg\t in %d calls to integrand\n",res,err, counter);

	res = integrate(&f_test2, a, b, abs, rel, &counter, &err);
	printf("int(4*sqrt(1-x*x))[0,1] = pi = \t %lg +/- %lg\t in %d calls to integrand\n\n",res,err, counter);




	double CC_test(double x, int* counter){
		(*counter)++;
		return 1.0/sqrt(x);
	}
	double CC_test2(double x, int* counter){
		(*counter)++;
		return log(x)/sqrt(x);
	}
	res = integrate(CC_test, 0, 1, abs, rel, &counter, &err);
	printf("without CC,\t int(1/sqrt(x))[0,1] = 2 =\t %lg +/- %lg\t in %d calls to integrand\n",res,err, counter);

	res = integrate_CC(CC_test, 0, 1, abs, rel, &counter, &err);
	printf("via CC,\t\t int(1/sqrt(x))[0,1] = 2 =\t %lg +/- %lg\t in %d calls to integrand\n",res,err, counter);

	res = integrate(CC_test2, 0, 1, abs, rel, &counter, &err);
	printf("without CC,\t int(ln(x)/sqrt(x))[0,1] = -4 =\t %lg +/- %lg\t in %d calls to integrand\n",res,err, counter);

	res = integrate_CC(CC_test2, 0, 1, abs, rel, &counter, &err);
	printf("via CC,\t\t int(ln(x)/sqrt(x))[0,1] = -4 =\t %lg +/- %lg\t in %d calls to integrand\n",res,err, counter);
	




	printf("\nDigits of pi via int(4*sqrt(1-x*x))[0,1] = PI:       for abs = rel = 1e-4 for all\n");
	printf("pi from memory = 3.141592653589793238462643383279502884197169399375105820\n");
	abs = 1e-4;
	rel = 1e-4;


	res = integrate(f_test2, 0, 1, abs, rel, &counter, &err);
	printf("without CC,\t\t  calls = %7d \t  => %.20lg \t |diff| = %.5lg\n", counter,res, fabs(3.141592653589793238462643383279502884197169399375105820-res));
	res = integrate_CC(f_test2, 0, 1, abs, rel, &counter, &err);
	printf("with CC,\t\t  calls = %7d \t  => %.20lg \t |diff| = %.5lg\n", counter,res, fabs(3.141592653589793238462643383279-res));


	//gsl_integration_workspace* integ = gsl_integration_workspace_alloc(1000);
	double result = 0;
	double gslerr = 0;
	size_t nev = 0;
	double gsl_f(double x, void* params){return 4*sqrt(1-x*x);}

	gsl_function F;
	F.function = &gsl_f;
	F.params = NULL;
	gsl_integration_qng(&F, 0, 1, abs, rel, &result, &gslerr, &nev);

	printf("with simple gsl,\t  calls = %7ld \t  => %.20lg \t |diff| = %.5lg\n", nev,result, fabs(3.141592653589793238462643383279-result));

	gsl_integration_workspace* integ = gsl_integration_workspace_alloc(1000);
	gsl_integration_qag(&F, 0, 1, abs, rel, 1000, 3, integ, &result, &gslerr);
	gsl_integration_workspace_free(integ);
	printf("with adaptive gsl,\t  calls = dunno \t  => %.20lg \t |diff| = %.5lg\n",result, fabs(3.141592653589793238462643383279-result));

	printf("\nEvery routine can fairly easily get within machine epsilon of M_PI\n");
	printf("With more time, it would be cool to make a plot with accuracy given to the integration along the x-axis,\n");
	printf("and deviation from PI along the y-axis, for the different routines to see how fast each of them reach PI\n");







	printf("\n\nEvauating infinite integrals:\n");

	abs = 1e-8;
	rel = 1e-8;

	double inf_test1(double x, int* counter){
		(*counter)++;
		return 1.0/(x*x);
	}
	double inf_test2(double x, int* counter){
		(*counter)++;
		return 2.0/(1+x*x);
	}
	double inf_test3(double x, int* counter){
		(*counter)++;
		return exp(-x*x);
	}
	res = integrate_inf(inf_test1, 1, INFINITY, abs, rel, &counter, &err);
	printf("Integral of 1/(x*x)[1,inf] = 1   \t\t= %lg +/- %lg \t in %d calls\n",res,err,counter);
	res = integrate_inf(inf_test2, 0, INFINITY, abs, rel, &counter, &err);
	printf("Integral of 2/(1+x*x)[0,inf] = pi \t\t= %lg +/- %lg \t in %d calls\n",res,err,counter);
	res = integrate_inf(inf_test3, -INFINITY, INFINITY, abs, rel, &counter, &err);
	printf("Integral of exp(-x*x)[-inf,inf] = sqrt(pi) \t= sqrt(%lg) +/- %lg \t in %d calls\n",res*res,err,counter);


	printf("\nThe numbers are correct, so the relevant comparison between my integration and GSL\n");
	printf("   would be the number of calls, but number of calls is not easily given by GSL integration\n");
	printf(" I would need to restructure the way my code works to allow gsl_function to also increment my counter\n");
	printf(" Since i have not done this, im not sure im eligible for the final point\n");



	return 0;
}