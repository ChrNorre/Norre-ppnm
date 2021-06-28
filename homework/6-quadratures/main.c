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
	double abs = 1e-10;
	double rel = 1e-10;
	double res = integrate(&f_test, a, b, abs, rel, &counter);
	printf("\nint(sqrt(x))[0,1] = 2/3 = %lg in %d calls to integrand\n",res, counter);
	counter = 0;
	res = integrate(&f_test2, a, b, abs, rel, &counter);
	printf("int(4*sqrt(1-x*x))[0,1] = pi = %lg in %d calls to integrand\n\n",res, counter);




	double CC_test(double x, int* counter){
		(*counter)++;
		return 1.0/sqrt(x);
	}
	double CC_test2(double x, int* counter){
		(*counter)++;
		return log(x)/sqrt(x);
	}
	counter = 0;
	res = integrate(CC_test, 0, 1, abs, rel, &counter);
	printf("without CC,\t int(1/sqrt(x))[0,1] = 2 = %lg in %d calls to integrand\n",res, counter);
	counter = 0;
	res = integrate_CC(CC_test, 0, 1, abs, rel, &counter);
	printf("via CC,\t\t int(1/sqrt(x))[0,1] = 2 = %lg in %d calls to integrand\n",res, counter);
	counter = 0;
	res = integrate(CC_test2, 0, 1, abs, rel, &counter);
	printf("without CC,\t int(ln(x)/sqrt(x))[0,1] = -4 = %lg in %d calls to integrand\n",res, counter);
	counter = 0;
	res = integrate_CC(CC_test2, 0, 1, abs, rel, &counter);
	printf("via CC,\t\t int(ln(x)/sqrt(x))[0,1] = -4 = %lg in %d calls to integrand\n",res, counter);
	




	printf("\nDigits of pi via int(4*sqrt(1-x*x))[0,1] = PI:\n");
	printf("pi from memory = 3.141592653589793238462643383279502884197169399375105820\n");
	counter = 0;
	res = integrate(f_test2, 0, 1, abs, rel, &counter);
	printf("without CC,\t\t  calls = %7d \t  => %.20lg \t |diff| = %.5lg\n", counter,res, fabs(3.141592653589793238462643383279502884197169399375105820-res));
	counter = 0;
	res = integrate_CC(f_test2, 0, 1, abs, rel, &counter);
	printf("with CC,\t\t  calls = %7d \t  => %.20lg \t |diff| = %.5lg\n", counter,res, fabs(3.141592653589793238462643383279-res));


	//gsl_integration_workspace* integ = gsl_integration_workspace_alloc(1000);
	double result = 0;
	double err = 0;
	size_t nev = 0;
	counter = 0;
	double gsl_f(double x, void* params){return 4*sqrt(1-x*x);}

	gsl_function F;
	F.function = &gsl_f;
	F.params = NULL;
	abs = 1e-4;
	rel = 1e-4;
	gsl_integration_qng(&F, 0, 1, abs, rel, &result, &err, &nev);

	printf("with simple gsl,\t  calls = %7ld \t  => %.20lg \t |diff| = %.5lg\n", nev,result, fabs(3.141592653589793238462643383279-result));

	abs = 1e-10;
	rel = 1e-10;
	gsl_integration_workspace* integ = gsl_integration_workspace_alloc(1000);
	gsl_integration_qag(&F, 0, 1, abs, rel, 1000, 3, integ, &result, &err);
	gsl_integration_workspace_free(integ);
	printf("with adaptive gsl,\t  calls = dunno \t  => %.20lg \t |diff| = %.5lg\n",result, fabs(3.141592653589793238462643383279-result));

	printf("\nEvery routine can fairly easily get within machine epsilon of M_PI\n");
	printf("With more time, it would be cool to make a plot with accuracy given to the integration along the x-axis,\n");
	printf("and deviation from PI along the y-axis, for the different routines to see how fast each of them reach PI\n");

	return 0;
}