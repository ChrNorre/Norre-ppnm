#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>


double et2(double x, void* params){
	double result = exp(-1*pow(x,2));
	return result;
}

double my_int_0z(double z, gsl_integration_workspace* w,int limit){
	double a=0,b=z,acc=1e-6,eps=1e-6,result,error;
	gsl_function F;
	F.function=&et2;
	gsl_integration_qags(&F, a, b, acc, eps, limit, w, &result, &error);
	result *= 2/(sqrt(M_PI));
	return  result;
}


int main(){
	int limit=999;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc (limit);

	for (double z=-3; z < 3; z += 1.0/8) {
		printf("%10g %10g %10g\n",z,my_int_0z(z,w,limit),erf(z));
	}
	
	gsl_integration_workspace_free(w);
return 0;
}