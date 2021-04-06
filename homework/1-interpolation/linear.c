int binsearch(int n, double* x, double z);

double linterp(int n, double* x, double* y, double z){
	int i = binsearch(n, x, z);
	return y[i] + (y[i+1]-y[i]) / (x[i+1]-x[i]) * (z-x[i]);
}


double linterp_integ(int n, double* x, double* y, double z){
	int i = binsearch(n, x, z);
	double res = 0;
	for (int j=0; j < i; j++){
		res += 0.5 * (x[j+1] - x[j]) * (y[j+1] + y[j]);
	}
	double zy = linterp(n, x, y, z);
	res += 0.5 * (z - x[i]) * (zy + y[i]);
	return res;
}