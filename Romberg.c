#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include "Romberg.h"

#define EPS 1.0e-8
#define JMAX 40
#define JMAXP (JMAX+1)
#define K 5

#define NRANSI
#define FREE_ARG char*
#define NR_END 1

#define FUNC(x) ((*func)(x))
double trapzd(double (*func)(double), double a, double b, int n)
{
    double x,tnm,sum,del;
    static double s;
    int it,j;

    if (n == 1) {
        return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
    } else {
        for (it=1,j=1;j<n-1;j++) it <<= 1;       //00000010 -> 00000100  (2^1->2^n)
        tnm=it;
        del=(b-a)/tnm;
        x=a+0.5*del;
        for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);      //integration in [a,b] for 2^n-2^1 bins
        s=0.5*(s+(b-a)*sum/tnm);
        return s;
    }
}
#undef FUNC

double qromb(double (*func)(double), double a, double b)
{
    double ss,dss;
    double s[JMAXP],h[JMAXP+1];
    int j;

    h[1]=1.0;
    for (j=1;j<=JMAX;j++) {
        s[j]=trapzd(func,a,b,j);                                   //return T_{2j}
        if (j >= K) {
            polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
            if (fabs(dss) <= EPS*fabs(ss)) return ss;
        }
        h[j+1]=0.25*h[j];
    }
    nrerror("Too many steps in routine qromb");
    return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K



void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
    int i,m,ns=1;
    double den,dif,dift,ho,hp,w;
    double *c,*d;

    dif=fabs(x-xa[1]);
    c=dvector(1,n);
    d=dvector(1,n);
    for (i=1;i<=n;i++) {
        if ( (dift=fabs(x-xa[i])) < dif) {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    *y=ya[ns--];
    for (m=1;m<n;m++) {
        for (i=1;i<=n-m;i++) {
            ho=xa[i]-x;
            hp=xa[i+m]-x;
            w=c[i+1]-d[i];
            if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
            den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }
    free_dvector(d,1,n);
    free_dvector(c,1,n);
}
#undef NRANSI


void nrerror(char error_text[])
    /* Numerical Recipes standard error handler */
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}

double *dvector(long nl, long nh)
    /* allocate a double vector with subscript range v[nl..nh] */
{
    double *v;

    v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
    if (!v) nrerror("allocation failure in dvector()");
    return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
    /* free a double vector allocated with dvector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

