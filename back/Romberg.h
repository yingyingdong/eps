
extern double qromb(double (*func)(double), double a, double b);
extern double test_func(double a);
extern double *dvector(long nl, long nh);
extern void free_dvector(double *v, long nl, long nh);
extern void nrerror(char error_text[]);
extern void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
extern double trapzd(double (*func)(double), double a, double b, int n);

