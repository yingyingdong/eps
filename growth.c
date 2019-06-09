#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Romberg.h"
#include "power.h"

double growth_int(double a)
{
    return pow(a / (Omega + (1 - Omega - OmegaLambda) * a + OmegaLambda * pow(a,-3*ww)), 1.5);
}

double growth(double a)
{
    double hubble_a;

    hubble_a = sqrt(Omega / (a * a * a) + (1 - Omega - OmegaLambda) / (a * a) + OmegaLambda*pow(a,-3*(1+ww)));

    //printf("hubble_a=%lf\n",hubble_a);
    return hubble_a * qromb(growth_int, 0, a);
}

double GrowthFactor(double astart, double aend)
{
    printf("%lf,%lf\n",growth(aend),growth(astart));
    return growth(aend) / growth(astart);
}



double growth_intdfy(double a)
{
    double omega_z,Ez2;
    Ez2 = Omega/a/a/a + (1-OmegaLambda-Omega)/a/a + OmegaLambda*pow(a,-3*(1+ww));
    omega_z = Omega*pow(a,-3)/Ez2;
    if(ww>-1) return 1./a* (pow(omega_z,0.55+0.05*(1+ww)) - 1);
    if(ww<-1) return 1./a* (pow(omega_z,0.55+0.02*(1+ww)) - 1);
}

double growth_dfy(double a)
{
    //printf("hubble_a=%lf\n",sqrt(omega/a/a/a + (1-omegalambda-omega)/a/a + omegalambda*pow(a,-3*(1+ww))));
    return exp(qromb(growth_intdfy,1e-3,a));
}

double GrowthFactor_dfy(double astart, double aend)
{
    printf("%lf,%lf\n",growth_dfy(aend),growth_dfy(astart));
    return (aend/astart)*growth_dfy(aend) / growth_dfy(astart);
}

double GF(double astart, double aend)
{
    if(ww==-1)
        return GrowthFactor(astart,aend);
    if(ww!=-1)
        return GrowthFactor_dfy(astart,aend);
}

double EZ_dfy(double a)
{
    return sqrt(Omega / (a * a * a) + (1 - Omega - OmegaLambda) / (a * a) + OmegaLambda*pow(a,-3*(1+ww)));
}

double hubble_dfy(double a)
{
    return Hubble*EZ_dfy(a);
}

double Omega_dfy(double a)
{
    return Omega*pow(a,-3)/EZ_dfy(a)/EZ_dfy(a);
}
