#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main(){
	
	//This is getting values from stdin
	double x;
	printf("x\t\tsin(x)\t\tcos(x)\n");
	
	int items;
	do {
		items = fscanf(stdin,"%lg",&x);
		printf("%lg\t\t%.4lg\t\t%.4lg\n",x,sin(x),cos(x));
	} while(items != EOF);
	
	
return 0;
}
