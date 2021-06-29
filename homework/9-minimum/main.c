#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include<stdio.h>
#include"minimum.h"


double Rosenbrock(gsl_vector* v){
	double x = gsl_vector_get(v,0);
	double y = gsl_vector_get(v,1);
	return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
}

double Himmelblau(gsl_vector* v){
	double x = gsl_vector_get(v,0);
	double y = gsl_vector_get(v,1);
	return pow(x*x+y-11,2) + pow(x+y*y-7,2);
}



//global variables to make them availiable in the deviation function without having them as parameters
gsl_vector* Evec;
gsl_vector* sigvec;
gsl_vector* delvec;


double Breit_Wigner(double E, double m, double Gamma, double A){
		double denom = (E-m)*(E-m) + Gamma*Gamma/4;
		return A/denom;
	}

double deviation(gsl_vector* v){
	double m = gsl_vector_get(v,0);
	double Gamma = gsl_vector_get(v,1);
	double A = gsl_vector_get(v,2);
	double res = 0;
	for (int i=0; i<30; i++){
		double F = Breit_Wigner(gsl_vector_get(Evec,i), m, Gamma, A);
		res += pow((F-gsl_vector_get(sigvec,i)),2) / pow(gsl_vector_get(delvec,i),2);
	}
	return res;
}



int main(){

	int counter = 0;
	double eps = 1e-8;
	gsl_vector* x = gsl_vector_alloc(2);
	gsl_vector* gx = gsl_vector_alloc(2);

	printf("\nRosenbrock minimum:\n(x0,y0) = (0,1)\n");
	gsl_vector_set(x,0,0);
	gsl_vector_set(x,1,1);
	qnewton(Rosenbrock, x, eps, &counter);
	printf("Minimum found at (x,y) = (%f,%f), in %d steps\n",gsl_vector_get(x,0),gsl_vector_get(x,1), counter);
	num_gradient(Rosenbrock, x, gx);
	printf("Function value is %f, Gradient here is (%f,%f)\n",Rosenbrock(x), gsl_vector_get(gx,0),gsl_vector_get(gx,1));

	printf("\nHimmelblau minimum:\n(x0,y0) = (1,1)\n");
	gsl_vector_set(x,0,1);
	gsl_vector_set(x,1,1);
	qnewton(Himmelblau, x, eps, &counter);
	printf("Minimum found at (x,y) = (%f,%f), in %d steps\n",gsl_vector_get(x,0),gsl_vector_get(x,1), counter);
	num_gradient(Himmelblau, x, gx);
	printf("Function value is %f, Gradient here is (%f,%f)\n",Himmelblau(x), gsl_vector_get(gx,0),gsl_vector_get(gx,1));





	//Higgs
	printf("\nHiggs data analysis:\n");
	int Ei;
	double sigi, deli;
	FILE* data = fopen("higgs.txt","r");
	Evec = gsl_vector_alloc(30);
	sigvec = gsl_vector_alloc(30);
	delvec = gsl_vector_alloc(30);


	//loading data and saving in three vectors
	int items;
	int i = 0;
	do {
		items = fscanf(data,"%d %lg %lg ",&Ei, &sigi, &deli);
		if (items == EOF) break;
		//printf("%d\t\t%lg\t\t%lg\n\n",Ei, sigi, deli);
		gsl_vector_set(Evec,i,Ei);
		gsl_vector_set(sigvec,i,sigi);
		gsl_vector_set(delvec,i,deli);
		i += 1;
	} while(1);
	fclose(data);


	gsl_vector* v = gsl_vector_alloc(3);	//v=(m,Gamma,A)
	gsl_vector_set(v,0,1);
	gsl_vector_set(v,1,1);
	gsl_vector_set(v,2,1);
	printf("With an initial guess of (m0,Gamma0,A0) = (1,1,1)\nThe deviation function gets minimized for:\n");
	qnewton(deviation, v, eps, &counter);
	printf("(m,Gamma,A) = (%lg, %lg, %lg)\n",gsl_vector_get(v,0), gsl_vector_get(v,1), gsl_vector_get(v,2));
	printf("The mass of the Higgs boson is m=%lg\n",gsl_vector_get(v,0));
	printf("Wikipedia says the real value is 125.3±0.6\n");
	printf("The width of the Higgs boson is Gamma=%lg\n",gsl_vector_get(v,1));





	gsl_vector_free(x);
	gsl_vector_free(gx);
	gsl_vector_free(v);
	gsl_vector_free(Evec);
	gsl_vector_free(sigvec);
	gsl_vector_free(delvec);

	return 0;
}