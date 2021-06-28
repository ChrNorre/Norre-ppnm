#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include<stdio.h>
#include"rootfinder.h"
#include"ode.h"


void testfun_1D(gsl_vector* x,gsl_vector* fx){
	// f(x) = 1-x^2,    should have roots in +-1
	gsl_vector_set(fx, 0, 1-pow(gsl_vector_get(x,0),2));
}

void testfun_2D(gsl_vector* x,gsl_vector* fx){
	// f1(x1,x2) = 2-x1^2-x2^2;
	// f2(x1,x2) = 2-x1^3-x2^3;
	// should have a solution at (x1,x2) = (1,1)
	gsl_vector_set(fx, 0, 2-pow(gsl_vector_get(x,0),2)-pow(gsl_vector_get(x,1),2));
	gsl_vector_set(fx, 1, 2-pow(gsl_vector_get(x,0),3)-pow(gsl_vector_get(x,1),3));
}

void Rosenbrock(gsl_vector* x,gsl_vector* fx){
	//f(x,y) = (1-x)^2 + 100*(y-x^2)^2;
	//0 = grad(f(x,y)) =
	//   2x - 400x(y-x^2) - 2 = 0  			f1(x1,x2)
	//          200y - 200x^2 = 0 			f2(x1,x2)
	// has solution (x,y) = (1,1)
	double xa = gsl_vector_get(x,0);
	double ya = gsl_vector_get(x,1);
	gsl_vector_set(fx, 0, 2*xa-400*xa*(ya-xa*xa)-2);
	gsl_vector_set(fx, 1, 200*ya-200*xa*xa);
}

double E = -1;		//global variable so i can change E in function without having it as a parameter
void Hydrogen(double r, gsl_vector* f, gsl_vector* dfdr){
		//radial Hydrogen atom
		//-(1/2)f'' -(1/r)f= E*f
		gsl_vector_set(dfdr,0,gsl_vector_get(f,1));       	       //y1' = y2
		gsl_vector_set(dfdr,1, -2*(E+1.0/r)*gsl_vector_get(f,0));   //y2' = -2(E + 1/r)*f
	}

int print_fun = 0;
void Hydrogen_system(gsl_vector* newE, gsl_vector* fres){
	//this is the function i give to rootfinder. I varies newE (1-dim) to find when fres=F(R=10) = 0, as desired

	gsl_vector* f0 = gsl_vector_alloc(2);
	gsl_vector* fR = gsl_vector_alloc(2);
	
	double a = 0.01, b = 10, h = 0.001;		//starting close to r=0
	double acc = 1e-8;
	double eps = 1e-8;

	gsl_vector_set(f0,0,a-a*a);	//F(r->0) = r-r^2 -> 0
	gsl_vector_set(f0,1,1-2*a);	//F'(r->0) = 1-2r -> 1


	E = gsl_vector_get(newE,0);
	if (print_fun){
		FILE* path = fopen("Hydrogen.txt","w");
		driver(&Hydrogen, a, f0, b, fR, h, acc, eps, path);
		fclose(path);
	}
	else {
		driver(&Hydrogen, a, f0, b, fR, h, acc, eps, NULL);
	}

	gsl_vector_set(fres,0,gsl_vector_get(fR,0));

	gsl_vector_free(f0);
	gsl_vector_free(fR);
}

int main(){

	printf("\nRootfinding\n");
	double eps = 1e-9;


	gsl_vector* x1 = gsl_vector_alloc(1);
	gsl_vector_set(x1,0,0.5);
	printf("\n1D testfunction: f(x)=1-x^2, \nx0 = 0.5\n");
	newton(testfun_1D, x1, eps);
	printf("Root found at x = %f, should be x=+-1\n",gsl_vector_get(x1,0));
	gsl_vector_set(x1,0,-0.5);
	printf("x0 = -0.5\n");
	newton(testfun_1D, x1, eps);
	printf("Root found at x = %f, should be x=+-1\n",gsl_vector_get(x1,0));



	gsl_vector* x2 = gsl_vector_alloc(2);
	gsl_vector_set(x2,0,2);
	gsl_vector_set(x2,1,3);
	printf("\n2D testfunction: f1(x,y)=2-x^2-y^2,  f2(x,y)=2-x^3-y^3\n");
	printf("x0 = 2, y0 = 3\n");
	newton(testfun_2D, x2, eps);
	printf("Root found at (x,y) = (%f,%f), should be (x,y) = (1,1)\n",gsl_vector_get(x2,0),gsl_vector_get(x2,1));



	gsl_vector_set(x2,0,2);
	gsl_vector_set(x2,1,3);
	printf("\nRosenbrock: f(x,y) = (1-x)^2 + 100*(y-x^2)^2\n");
	printf("Finding extremum by rootfinding on its gradient\n");
	printf("f1(x,y)=2x - 400x(y-x^2) - 2 = 0\n");
	printf("f2(x,y)=200y - 200x^2 = 0\n");
	printf("x0 = 2, y0 = 3\n");
	newton(Rosenbrock, x2, eps);
	printf("Root found at (x,y) = (%f,%f), should be (x,y) = (1,1)\n",gsl_vector_get(x2,0),gsl_vector_get(x2,1));





	gsl_vector* Enew = gsl_vector_alloc(1);
	gsl_vector_set(Enew,0,-0.7);		
	// Strange things happen when starting guess is greater (more positive) than -0.5
	// Seeing as i was tasked to find the lowest energy, it would make the most sense to have a starting guess
	// of 0, and having it become more negative until a solution was found.
	newton(Hydrogen_system, Enew, eps);
	printf("\nSolving for minimum energy of radial part of hydrogen atom\n");
	printf("E0 = -0.7,  E_root = %f, should be E_actual = -0.5\n",gsl_vector_get(Enew,0));

	print_fun = 1;
	Hydrogen_system(Enew, Enew);



	gsl_vector_free(x1);
	gsl_vector_free(x2);
	gsl_vector_free(Enew);
	return 0;
}