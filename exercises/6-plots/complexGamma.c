#include<math.h>
#include<stdio.h>
#include<complex.h>

double complex Gamma(double complex z){
	///single precision gamma function (Gergo Nemes, from Wikipedia)
	if(creal(z)<0)return M_PI/csin(M_PI*z)/Gamma(1-z);
	if(creal(z)<9)return Gamma(z+1)/z;
	complex lnGamma=z*clog(z+1/(12*z-1/(z*10)))-z+clog(2*M_PI/z)/2;
	return cexp(lnGamma);
}


int main(){
	double width = 4.2;
	FILE* xyz_out = fopen("xyz_out.txt","w");


	for (double i = -width; i < width; i += 1.0/32){
		for(double j = -width; j < width; j += 1.0/32){
			fprintf(xyz_out, "%8g %8g %8g\n ", i, j, fmin(6,cabs(Gamma(i + j*I))));
			
		}
	}

	fclose(xyz_out);


return 0;
}