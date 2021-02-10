#include<stdio.h>
//#include<complex.h>
//#include<math.h>
#include<tgmath.h>
#define pi M_PI
#define e M_E


int main(){
	printf("\nGamma(5) = %g\n",tgamma(5));
	printf("Bessel_1(0.5) = %g\n",j1(0.5));
	complex z1 = csqrt(-2);
	printf("sqrt(-2) = %g + I%g\n",creal(z1),cimag(z1));
	complex z2 = cexp(I*pi);
	printf("exp(i*pi) = %g + I%g\n",creal(z2),cimag(z2));
	complex z3 = cexp(I);
	printf("exp(i) = %g + I%g\n",creal(z3),cimag(z3));
	complex z4 = cpow(e,I);
	printf("exp(i) = cpow(e,i) = %g + I%g\n",creal(z4),cimag(z4));
	complex z5 = cpow(I,e);	
	printf("cpow(i,e) = %g + I%g\n",creal(z5),cimag(z5));
	complex z6 = cpow(I,I);	
	printf("cpow(i,i) = %g + I%g\n",creal(z6),cimag(z6));

	float x_f = 1.f/9;
	double x_d = 1./9;
	long double x_ld = 1.L/9;
	printf("\nfloat decimals:\t \t %.25g\n", x_f);
	printf("double decimals:\t %.25lg\n", x_d);
	printf("long double decimals:\t %.25Lg\n", x_ld);
	


return 0;
}
