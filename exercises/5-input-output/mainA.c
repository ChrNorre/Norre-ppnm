#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main(int argc, char** argv){
	
	//This is getting values from arguments to the executable
	double x;
	printf("x\t\tsin(x)\t\tcos(x)\n");
	if(argc>1) {
		for (int i= 1; i<argc; i++) {
			x = atof(argv[i]);
			printf("%lg\t\t%.4lg\t\t%.4lg\n",x,sin(x),cos(x));
		} 
	}

return 0;
}