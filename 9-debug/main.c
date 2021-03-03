#include"stdio.h"
//#include"gsl_matrix.h" //usualy in gsl/gsl_something
#include"gsl/gsl_matrix.h"

int print_half_00(gsl_matrix* m)
{
	//double half = 1/2; 				//integer division
	double half = 1.0/2;
	//int *status = printf( "half m_{00} = %i\n", gsl_matrix_get(&m,0,0)*half ); 
	//m is already a pointer, so & is not needed. %i is a int placeholder, but matrix_get and half are doubles (use %g). return type of printf is just int
	int status = printf( "half m_{00} = %g\n", gsl_matrix_get(m,0,0)*half );
	//gsl_matrix_free(m); likely better to free in main
	return status;
}

int main(void)
{
	//gsl_matrix m = gsl_matrix_alloc(0,0);		//alloc returns a pointer, size should be 1,1 to be able to index at 0,0
	gsl_matrix* m = gsl_matrix_alloc(1,1);
	//gsl_matrix_set(&m,0,0,66);		//m is already a pointer
	gsl_matrix_set(m,0,0,66);
	printf("half m_{00} (should be 33):\n");
	//int *status = print_half_00(&m); //m is already a pointer
	int status = print_half_00(m);
	//if(status>0)  		//printf failure returns a negative number
	if(status<0)
		//printf("status=%g : SOMETHING WENT TERRIBLY WRONG (status>0)\n",*status); //status is int, so use %i
		printf("status=%i : SOMETHING WENT TERRIBLY WRONG (status<0)\n",status);
	else
		//printf("status=%g : everything went just fine (status=0)\n",*status); //same as above
		printf("status=%i : everything went just fine (status=0)\n",status);
	//gsl_matrix_free(&m); //m is already a pointer
	gsl_matrix_free(m);
return 0;
}