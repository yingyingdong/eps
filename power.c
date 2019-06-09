#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Romberg.h"
#include "power.h"

static double r_tophat;
static double r_sharpk;
static double Norm;
static double R8;
double res;

double PowerSpec(double k)
{
    double power = 0, alpha, Tf;

    power = PowerSpec_EH(k);
    power *= pow(k, PrimordialIndex - 1.0);
    return power;
}

double PowerSpec_DM_2ndSpecies(double k)
{
    /* at the moment, we simply call the Eistenstein & Hu spectrum
     * for the second DM species, but this could be replaced with
     *  something more physical, say for neutrinos
     * */

    double power; 
    power = Norm * k * pow(tk_eh(k), 2);
    return power;
}


double PowerSpec_EH(double k)   /* Eisenstein & Hu */
{
    return Norm * k * pow(tk_eh(k), 2);
}


double tk_eh(double k)      /* from Martin White */         //EH transfer function
{
    double q, theta, ommh2, a, s, gamma, L0, C0;
    double tmp;
    double omegam, ombh2, hubble;

    /* other input parameters */
    hubble = HubbleParam;

    omegam = Omega;
    ombh2 = OmegaBaryon * HubbleParam * HubbleParam;

    if(OmegaBaryon == 0)
        ombh2 = 0.04 * HubbleParam * HubbleParam;

    k *= (3.085678e24 / UnitLength_in_cm);    /* convert to h/Mpc */

    theta = 2.728 / 2.7;
    ommh2 = omegam * hubble * hubble;
    s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
    a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
        + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
    gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
    gamma *= omegam * hubble;
    q = k * theta * theta / gamma;
    L0 = log(2. * exp(1.) + 1.8 * q);
    C0 = 14.2 + 731. / (1. + 62.5 * q);
    tmp = L0 / (L0 + C0 * q * q);
    return (tmp);
}

double TopHatSigma2(double R)
{
    r_tophat = R;
    r_sharpk = R;
//  printf("TopHatSigma2:R=%12g,Norm=%12g\n",R,Norm);
    if(tophat==1) return qromb(sigma2_int_tophat, 0, 500.0 * 1 / R);   /* note: 500/R is here chosen as integration boundary (infinity) */
    else return qromb(sigma2_int_sharpk, 0, 1./R);
}

double sigma2_int_tophat(double k) //Top Hat (TH) window function (in real space)
{
    double kr, kr3, kr2, w, x;

    kr  = r_tophat * k;
    kr2 = kr * kr;
    kr3 = kr2 * kr;

    if(kr < 1e-8)
        return 0;

    w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);        //window function
    x = 4 * PI * k * k * w * w * PowerSpec(k);      //P207, (4.271)

    return x;
}

double sigma2_int_sharpk(double k) //sharp k-space (SK) window function (in k space)
{
    double kr,kr3,kr2,w,x;
    
    if(k<1./r_sharpk) w = 1;
    else              w = 0;

    x = 4 * PI * k * k *w * w * PowerSpec(k);
    return x;
}

void initialize_powerspectrum(void)
{
    Norm = 1.0;
    R8 = 8 * (3.085678e24 / UnitLength_in_cm);    /* 8 Mpc/h */
    res = TopHatSigma2(R8);

    printf("\nNormalization of spectrum in file:  Sigma8 = %g\n", sqrt(res));
    Norm = Sigma8 * Sigma8 / res;
    printf("\nNormalization of spectrum in file:  Sigma8 = %g,Norm=%lg\n", sqrt(res),Norm);
}
