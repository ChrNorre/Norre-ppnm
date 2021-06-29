double integrate(double (*f)(double, int*), double a, double b, double abs, double rel, int* counter, double* err);
double integrate24(double (*f)(double, int*), double a, double b, double abs, double rel, double p26, double p46, int* counter, double* err);
double integrate_CC(double (*f)(double, int*), double a, double b, double abs, double rel, int* counter, double* err);
double integrate24_CC(double (*CC)(double, int*), double a, double b, double abs, double rel, double p26, double p46, int* counter, double* err);
double integrate_inf(double (*f)(double, int*), double a, double b, double abs, double rel, int* counter, double* err);