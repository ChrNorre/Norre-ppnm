#include<stdio.h>
#include<limits.h>
#include<float.h>

//i only have to specify the function header. the implementation is handled in the linking.
int equal(double a, double b, double tau, double epsilon);


void int_max(){
	// I use gcc, and my i continues to grow larger than INT_MAX, (does not overflow to negative numbers)
	for(int i= INT_MAX-10; i+1>i; i++){
		printf("INT_MAX - i = %d\n",INT_MAX-i);
	}

	int i = INT_MAX-10;
	while(i+1>i){
		printf("INT_MAX - i = %d\n",INT_MAX-i);
		i++;
	}

	i = INT_MAX-10;
	do {
		printf("INT_MAX - i = %d\n",INT_MAX-i);
		i++;
	} while(i+1>i);
}

void int_min(){
	//same problem as INT_MAX, 
	
	int i = INT_MIN+10;
	while(i-1<i){
		i--;
	}
	printf("min int is = %d\n",i);

	for(int i=INT_MIN+10; i-1<i; i--) {
		printf("INT_MIN - i = %d\n",INT_MIN-i);
	}
	
	i = INT_MIN+10;
	do {i--;} while(i-1<i);
	printf("min int is = %d\n",i);

}

void find_epsilon(){
	{double x = 1;
	while (1+x != 1) {x /= 2;}
	printf("\nDouble epsilon = %.14g\n",x*2);
	printf("DBL_EPSILON = %.14g\n",DBL_EPSILON);}

	
	{float x = 1;
	while (1+x != 1) {x /= 2;}
	printf("\nFloat epsilon = %.14f\n",x*2);
	printf("FLT_EPSILON = %.14f\n",FLT_EPSILON);}

	{long double x = 1;
	while (1+x != 1) {x /= 2;}
	printf("\nLong Double epsilon = %.14Lg\n",x*2);
	printf("LDBL_EPSILON = %.14Lg\n",LDBL_EPSILON);}

	//for loops
	{double x;
	for(x = 1; 1+x != 1; x /= 2) {}	
	x *=2;}
	{float x;
	for(x = 1; 1+x != 1; x /= 2) {}
	x *=2;}
	{long double x;
	for(x = 1; 1+x != 1; x /= 2) {}
	x *=2;}

	//do-while
	{double x=1;
	do {x /= 2;} while(1+x != 1);
	x *= 2;}
	{float x=1;
	do {x /= 2;} while(1+x != 1);
	x *= 2;}
	{long double x=1;
	do {x /= 2;} while(1+x != 1);
	x*= 2;}

}

void float_sum(){
	float sum_up_float = 0;
	float sum_down_float = 0;
	for(int i = 1; i<INT_MAX/3; i++) {
		sum_up_float += 1.0f/i;
		sum_down_float += 1.0f/(INT_MAX/3 - i);
	}
	printf("\nsum_up_float = %.8f\n",sum_up_float);
	printf("sum_down_float = %.8f\n",sum_down_float);
	printf("up-down = %.8f\n",sum_up_float-sum_down_float);
	/*
	 * sum_down_float is larger (by 3.something)
	 * The only explaination i can come up with is that the numbers are stored in scientific notation,
	 *  such that if we add two numbers of the same magnitude, there is greater precision.
	 *  sum_up_float is (something increasing) + (something decreasing), this will eventually 
	 *  have a difference in magnitude so large, that the smaller gets beyond FLT_EPSILON.
	 *  sum_down_float is (something increasing) + (something increasing), these will always be of approximately
	 *  the same magnitude, so there should be greater precision.
	 *  The difference does not depend on max, (at least no difference between /50 and /2)
	 */

	double sum_up_double = 0;
	double sum_down_double = 0;
	for(int i = 1; i<INT_MAX/2; i++){
		sum_up_double += 1.0f/i;
		sum_down_double += 1.0f/(INT_MAX/2 - i);

	}
	printf("\nsum_up_double = %.14g\n",sum_up_double);
	printf("sum_down_double = %.14g\n",sum_down_double);
	printf("up-down = %.14g\n",sum_up_double-sum_down_double);
	//double almost fixes it (diff = 1e-9)
	//This means that 1/INT_MAX and DBL_EPSILON is of approximately the same magnitude 

}


int main(){
	
	//int_max();
	//int_min();
	find_epsilon();
	float_sum();
	

	double a = 1.01;
	double b = 2.11;
	double tau = 1;
	double epsilon = 1;
	printf("\nequaltest = %d\n",equal(a, b, tau, epsilon));

	
return 0;
}
