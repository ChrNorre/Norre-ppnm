#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<math.h>


void driver(
	void (*f)(double,gsl_vector*,gsl_vector*), 
	double a, gsl_vector* ya, 
	double b, gsl_vector* yb, 
	double h, 
	double acc, double eps,
	FILE* path);


void HO(double t, gsl_vector* y, gsl_vector* dydt){
	//harmonic oscilaltor
	//u'' = -u
	gsl_vector_set(dydt,0,gsl_vector_get(y,1));        //y1' = y2
	gsl_vector_set(dydt,1,-1.0*gsl_vector_get(y,0));   //y2' = f(x,y1,y2) = -y1
}

void HO_system(){
	gsl_vector* ya = gsl_vector_alloc(2);
	gsl_vector* yb = gsl_vector_alloc(2);

	gsl_vector_set(ya,0,1);
	gsl_vector_set(ya,1,0);

	double a = 0, b = 5, h = 0.1;
	double acc = 1e-10;
	double eps = 1e-10;

	driver(&HO, a, ya, b, yb, h, acc, eps, NULL);

	printf("\nHarmonic oscillator from t=0 to t=5:\n");
	printf("x(0)=1, x(5)=%.7g\n",gsl_vector_get(yb,0));
	printf("v(0)=0, v(5)=%.7g\n",gsl_vector_get(yb,1));

	gsl_vector_free(ya);
	gsl_vector_free(yb);
}

void SIR(double t, gsl_vector* y, gsl_vector* dydt){
	//SIR model
	int N = 5.8e6;
	double Tc = 1.0;
	double Tr = 14.0;

	gsl_vector_set(dydt, 0, -gsl_vector_get(y,1)*gsl_vector_get(y,0)/(N*Tc)); //dS = -IS/NTc
	gsl_vector_set(dydt, 1, gsl_vector_get(y,1)*(gsl_vector_get(y,0)/(N*Tc) -1/Tr)); //dI = I(S/NTc - 1/Tr)
	gsl_vector_set(dydt, 2, gsl_vector_get(y,1)/Tr);  //dR = I/Tr
}

void SIR_system(){
	gsl_vector* ya = gsl_vector_alloc(3);
	gsl_vector* yb = gsl_vector_alloc(3);

	gsl_vector_set(ya,0,5.8e6);
	gsl_vector_set(ya,1,10);
	gsl_vector_set(ya,2,0);

	double a = 0, b = 20, h = 0.1;
	double acc = 1e-6;
	double eps = 1e-6;

	FILE* path = fopen("SIR_path.txt","w");	
	driver(&SIR, a, ya, b, yb, h, acc, eps, path);
	fclose(path);

	printf("\nSIR system for denmark with 10 initial infected, and Tc = 1 (days), Tr = 14 (days):\n");
	printf("After 20 days:\n");
	printf("S = %8d\n",(int) gsl_vector_get(yb,0));
	printf("I = %8d\n",(int) gsl_vector_get(yb,1));
	printf("R = %8d\n",(int) gsl_vector_get(yb,2));

	gsl_vector_free(ya);
	gsl_vector_free(yb);
}

int main(){

	HO_system();
	SIR_system();

	
	return 0;
}