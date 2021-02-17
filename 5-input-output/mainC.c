#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main(int argc, char** argv){
	
	//This is getting input-dir and output-dir as arguments
	//otherwise, similar to mainB.c
	
	FILE* infile = fopen(argv[1],"r");
	FILE* outfile = fopen(argv[2],"w");

	double x;
	fprintf(outfile,"x\t\tsin(x)\t\tcos(x)\n");
	int items;
	do {
		items = fscanf(infile,"%lg",&x);
		fprintf(outfile,"%lg\t\t%.4lg\t\t%.4lg\n",x,sin(x),cos(x));
	} while(items != EOF);

return 0;
}