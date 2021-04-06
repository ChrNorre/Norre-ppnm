#include <stdio.h>
double f(int i){
	return 1.0/i;
}

int main(){
	int n = (int) 1e7;
	printf("n=%g\n",(double) n);
	double s = 0;
	for(int i = 1; i < n; i++) {
		s += f(i);
	}
	printf("s=%g\n",s);
}

