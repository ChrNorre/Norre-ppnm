#include"ann.h"
#include"minimum.h"
#include<math.h>
#include<stdio.h>
#define get gsl_vector_get
#define set gsl_vector_set



ann* ann_alloc(int n,double(*f)(double)){
	gsl_vector* p = gsl_vector_alloc(3*n);
	ann* network = malloc(sizeof(ann));
	network->n = n;
	network->f = f;
	network->params = p;
	return network;
}

void   ann_free(ann* network){
	gsl_vector_free(network->params);
	free(network);
}


double ann_response(ann* network, double x){

	int n = network->n;
	double(*f)(double) = network->f;
	gsl_vector* params = network->params;

	double ai, bi, wi, res = 0;
	for (int i=0; i<n; i++){
		ai = get(params, 3*i+0);
		bi = get(params, 3*i+1);
		wi = get(params, 3*i+2);
		res += (*f)((x-ai)/bi)*wi;
	}
	return res;
}



//field variables is not an elegant solution, especially outside of main.c, but i ran out of 
// other options
gsl_vector* x_field;
gsl_vector* y_field;
ann* network_field;

double C(gsl_vector* p){
	double res = 0;
	for (int j=0; j<x_field->size; j++){
		double Fp_x = ann_response(network_field, get(x_field,j));
		res += (Fp_x - get(y_field,j)) * (Fp_x - get(y_field,j));
	}
	return res;
}

void ann_train(ann* network,gsl_vector* xs,gsl_vector* ys){
	//the aim is to call this function and get a trained network back
	//this function needs to implement minimization of C(p)

	x_field = xs;
	y_field = ys;
	network_field = network;

	//I would like to define the C(p) function here, but then i get segmentation fault when passing it to qnewton
	// (specifically, the first time i try to use C(p))

	int counter = 0;
	double eps = 0.001;

	//fprintf(stderr,"C(p) before minimize = %lg\n",C(network->params));

	qnewton(C, network->params, eps, &counter);

	//fprintf(stderr,"C(p) after minimize = %lg\n",C(network->params));
	
	//fprintf(stderr,"done minimizing, steps=%d\n",counter);
	
	//for (int i=0; i < network->params->size; i++){
	//	fprintf(stderr,"p(%d)=%lg\n",i,get(network->params,i));
	//}

}


double ann_integ(ann* network, double (*F)(double), double z){
	// F is integral of activation function
	double res = 0;
	for(int i=0; i < network->n; i++){
		double ai = get(network->params,3*i+0);
		double bi = get(network->params,3*i+1);
		double wi = get(network->params,3*i+2);
		res += wi*bi*F((z-ai)/bi);		//from chain rule
	}
	return res;
}

double ann_deriv(ann* network, double (*df)(double), double z){
	// df is derivative of activation function
	double res = 0;
	for(int i=0; i < network->n; i++){
		double ai = get(network->params,3*i+0);
		double bi = get(network->params,3*i+1);
		double wi = get(network->params,3*i+2);
		res += wi*df((z-ai)/bi)/bi;		//from chain rule
	}
	return res;
}
