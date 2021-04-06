#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>

double Erf(double x){
	/// single precision error function (Abramowitz and Stegun, from Wikipedia)
	if(x<0) return -Erf(-x);
	double a[]={0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429};
	double t=1/(1+0.3275911*x);
	double sum=t*(a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*a[4]))));/* the right thing */
	return 1-sum*exp(-x*x);
}

int main(){
	double xmin = -2;
	double xmax = 2;
	for (double x=xmin; x<xmax; x += 1.0/8){
		printf("%8g %8g %8g %8g\n",x,erf(x),gsl_sf_erf(x), Erf(x));
	}
return 0;	
}


