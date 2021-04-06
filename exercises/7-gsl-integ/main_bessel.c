#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>


double cossin(double tau, void* params){
	/* This is my solution to multiple parameters (doesnt work for different types)
	I put my parameters (order and x) into a double array, then i take the pointer to that array,
	and casts it as a void pointer. Then, inside the function, i cast the void pointer back to a double pointer
	and then i add a number to it (since double pointers and double array pointers are the same)
	This has the effect of jumping one 'double' forward in memory, and if the original structure
	was an array, the value will make sense. I then read what is at that location in memory.
	The correct thing would probaby be to create a structure for the params.
	*/
	double n = *(((double *)params)+0);
	double x = *(((double *)params)+1);
	double result = cos(n*tau-x*sin(tau));
	return result;
}

double bessel_int(double n, double x, gsl_integration_workspace* w,int limit){
	double a=0,b=M_PI,acc=1e-6,eps=1e-6,result,error;
	gsl_function F;
	F.function=&cossin;
	double para[2] = {n,x};
	F.params = (void*) &para;
	gsl_integration_qags(&F, a, b, acc, eps, limit, w, &result, &error);
	result *= 1/M_PI;
	return  result;
}


int main(){
	int limit=999;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc (limit);

	FILE* bessel = fopen("bessel.txt","w");
	FILE* bessel_ref = fopen("bessel_ref.txt","w");

	for (double x=1.0/8; x < 15; x += 1.0/8) {
		fprintf(bessel,"%10g %10g %10g %10g\n",x,bessel_int(0.0,x,w,limit),bessel_int(1.0,x,w,limit),bessel_int(2.0,x,w,limit));
	}

	for (double x=1.0/2; x<15; x+= 1.0/2){
		fprintf(bessel_ref, "%10g %10g %10g %10g\n",x,jn(0,x),jn(1,x),jn(2,x));
	}

	fclose(bessel);
	fclose(bessel_ref);


	gsl_integration_workspace_free(w);
return 0;
}