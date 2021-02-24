#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>


double lnsqrt(double x, void* params){
	double nume = log(x);
	double deno = sqrt(x);
	return nume/deno;
}

double my_int_01(gsl_integration_workspace* w,int limit){
	double a=0,b=1,acc=1e-6,eps=1e-6,result,error;
	gsl_function F;
	F.function=&lnsqrt;
	gsl_integration_qags(&F, a, b, acc, eps, limit, w, &result, &error);
	return result;
}


int main(){
	int limit=999;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc (limit);


	printf("int of ln/sqrt = %g\n",my_int_01(w,limit));


	gsl_integration_workspace_free(w);
return 0;
}