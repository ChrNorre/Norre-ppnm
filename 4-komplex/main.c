#include"komplex.h"
#include<stdio.h>

int main(){
	komplex a = {1,2}, b = {3,4};
	komplex r;
	komplex R;


	printf("testing komplex_add...\n");
	r = komplex_add(a,b);
	R = komplex_new(4,6);
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a+b should   = ", R);
	komplex_print("a+b actually = ", r);

	printf("\ntesting komplex_sub...\n");
	r = komplex_sub(a,b);
	R = komplex_new(-2,-2);;
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a-b should   = ", R);
	komplex_print("a-b actually = ", r);

	printf("\ntesting komplex_conjugate...\n");
	r = komplex_conjugate(a);
	R = komplex_new(1,-2);;
	komplex_print("a=",a);
	komplex_print("a* should   = ", R);
	komplex_print("a* actually = ", r);

	printf("\ntesting komplex_abs...\n");
	double r_d = komplex_abs(b);
	double R_d = 5;
	komplex_print("b=",b);
	printf("abs(b) should   = %g\n", R_d);
	printf("abs(b) actually = %g\n", r_d);

	printf("\ntesting komplex_mul...\n");
	r = komplex_mul(a,b);
	R = komplex_new(-5,10);;
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a*b should   = ", R);
	komplex_print("a*b actually = ", r);

	printf("\ntesting komplex_div...\n");
	r = komplex_div(b,a);
	R = komplex_new(2.2,-0.4);;
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("b/a should   = ", R);
	komplex_print("b/a actually = ", r);
}
