#include<math.h>
#include<stdio.h>

double ex(double x){
	if(x<0)return 1/ex(-x);
	if(x>1./8)return pow(ex(x/2),2);
	return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}

double ex_neg(double x){
	if(x>1./8)return pow(ex(x/2),2);
	return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}

double ex_nohalf(double x){
	if(x<0)return 1/ex(-x);
	return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}

double ex_longsum(double x){
	if(x<0)return 1/ex(-x);
	if(x>1./8)return pow(ex(x/2),2);
	double res = 1;
	double run = 1;
	for (int i=1; i <=10; i++){
		run *= x/i;
		res += run;
	}
	return res;
}

double ex_naive(double x){
	double res = 1;
	double run = 1;
	for (int i=1; i <=10; i++){
		run *= x/i;
		res += run;
	}
	return res;
}


int main(){
	
	double x;
	for (int i = 0; i < 1000; i++){
		x = -3 + i*(13./1000);
		printf("%10g\t%10g\t%10g\t%10g\t%10g\t%10g\t%10g\t%10g\n",x, ex(x), exp(x),exp(x)-ex(x),fabs(ex(x)-ex_neg(x)),fabs(ex(x)-ex_nohalf(x)),fabs(ex(x)-ex_longsum(x)),fabs(ex(x)-ex_naive(x)));
	}




return 0;
}