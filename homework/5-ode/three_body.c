#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include"ode.h"

void print_vectors(gsl_vector* a, gsl_vector* b){
	for (int i=0; i<a->size; i++){
		printf("%10g %10g\n",gsl_vector_get(a,i),gsl_vector_get(b,i));
	}
}


void three_body(double t, gsl_vector* y, gsl_vector* dydt){
	//Threebody model
	gsl_vector_set(dydt, 0, gsl_vector_get(y,6));	//x1' = vx1
	gsl_vector_set(dydt, 1, gsl_vector_get(y,7));	//y1' = vy1
	gsl_vector_set(dydt, 2, gsl_vector_get(y,8));	//x2' = vx2
	gsl_vector_set(dydt, 3, gsl_vector_get(y,9));	//y2' = vy2
	gsl_vector_set(dydt, 4, gsl_vector_get(y,10));	//x3' = vx3
	gsl_vector_set(dydt, 5, gsl_vector_get(y,11));	//y3' = vy3

	double x1 = gsl_vector_get(y,0);
	double y1 = gsl_vector_get(y,1);
	double x2 = gsl_vector_get(y,2);
	double y2 = gsl_vector_get(y,3);
	double x3 = gsl_vector_get(y,4);
	double y3 = gsl_vector_get(y,5);

	double r12 = sqrt(pow((x1-x2),2) + pow((y1-y2),2));
	double r23 = sqrt(pow((x2-x3),2) + pow((y2-y3),2));
	double r31 = sqrt(pow((x3-x1),2) + pow((y3-y1),2));

	gsl_vector_set(dydt, 6, (x2-x1)/pow(r12,3) + (x3-x1)/pow(r31,3)); 
	gsl_vector_set(dydt, 7, (y2-y1)/pow(r12,3) + (y3-y1)/pow(r31,3)); 
	gsl_vector_set(dydt, 8, (x3-x2)/pow(r23,3) + (x1-x2)/pow(r12,3)); 
	gsl_vector_set(dydt, 9, (y3-y2)/pow(r23,3) + (y1-y2)/pow(r12,3)); 
	gsl_vector_set(dydt, 10, (x2-x3)/pow(r23,3) + (x1-x3)/pow(r31,3)); 
	gsl_vector_set(dydt, 11, (y2-y3)/pow(r23,3) + (y1-y3)/pow(r31,3)); 
}

int main(){

	gsl_vector* veca = gsl_vector_alloc(12);
	gsl_vector* vecb = gsl_vector_alloc(12);


	//from footnote on wikipedia
	double xa = -0.97000436;
	double ya = 0.24308753;
	double vxa = 0.4662036850;
	double vya = 0.4323657300;
	double vxa2 = -0.93240737;
	double vya2 = -0.86473146;

	gsl_vector_set(veca,0,xa);
	gsl_vector_set(veca,1,ya);
	gsl_vector_set(veca,2,0);
	gsl_vector_set(veca,3,0);
	gsl_vector_set(veca,4,-xa);
	gsl_vector_set(veca,5,-ya);
	gsl_vector_set(veca,6,vxa);
	gsl_vector_set(veca,7,vya);
	gsl_vector_set(veca,8,vxa2);
	gsl_vector_set(veca,9,vya2);
	gsl_vector_set(veca,10,vxa);
	gsl_vector_set(veca,11,vya);

	double a = 0, b = 6.3259, h = 0.1;
	double acc = 1e-8;
	double eps = 1e-8;

	//FILE* path = fopen("SIR_path.txt","w");	
	driver(&three_body, a, veca, b, vecb, h, acc, eps, NULL);
	//fclose(path);
	printf("\nSimulating figure-8 three body system for 1 period\n");
	printf("t=0	     t=6.3259\n");
	print_vectors(veca, vecb);

	gsl_vector_free(veca);
	gsl_vector_free(vecb);


	
	return 0;
}