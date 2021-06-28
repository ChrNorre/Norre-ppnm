#include"ann.h"
#include<math.h>
#define get gsl_vector_get
#define set gsl_vector_set


int main(){

	int n = 8;

	//activation function, gaussian wavelet
	double f(double x){
		return x*exp(-x*x);
	}
	double F(double x){
		return -0.5*exp(-x*x);
	}
	double df(double x){
		return (1-2*x*x)*exp(-x*x);
	}

	ann* network = ann_alloc(n, f);

	//initializing parameters
	for (int i=0; i<n; i++){
		set(network->params,3*i+0, -1+i*2/(n-1));	//initial offset
		set(network->params,3*i+1, 0.1);	//initial scaling
		set(network->params,3*i+2, 0.5);	//initial weight
	}



	double data_function(double x){
		return sin(5*x-1)*exp(-x*x);
	}
	//double data_function_integ(double x){
	//	return 0; //Something horrible i dont want to deal with
	//}
	double data_function_deriv(double x){
		return (5*cos(5*x-1) - 2*x*sin(5*x-1))*exp(-x*x);
	}

	//constructing datapoints to interpolate
	int N = 40;
	gsl_vector* xs = gsl_vector_alloc(N);	
	gsl_vector* ys = gsl_vector_alloc(N);
	double x = -1;
	double delx = 2.0/(N-1);
	for(int i=0;i<N;i++){
		set(xs,i,x);
		set(ys,i,data_function(x));
		x += delx;
		printf("%lg %lg\n",get(xs,i),get(ys,i));
	}
	printf("\n\n");


	//training network in parameters to make best interpolation
	ann_train(network, xs, ys);



	//making data for plotting
	double z = get(xs,0);
	while (z < get(xs,N-1)){
		double yz = ann_response(network, z);
		double Iz = ann_integ(network, F, z);
		double dz = ann_deriv(network, df, z);
		double point_dz = data_function_deriv(z);
		printf("%lg %lg %lg %lg %lg\n",z,yz,Iz,dz,point_dz);
		z += 0.02;
	}
	printf("\n\n");

	gsl_vector_free(xs);
	gsl_vector_free(ys);
	ann_free(network);


	return 0;
}