#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Romberg.h"
#include "power.h"
#include "growth.h"

void set_units(void)        /* ... set some units */
{
    UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;

    G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
    Hubble = HUBBLE * UnitTime_in_s;                                               //~0.9972, unit in km/s h/kpc
    printf("Hublle=%lf,UnitTime_in_s=%lf,G=%lf\n",Hubble,UnitTime_in_s,G);
}

int prepare()
{
    set_units();
    initialize_powerspectrum(); //normalize according to sigma8
    return 0;
}

typedef struct
{
    int n;
    double *z1;
    double *DZ;
    double *M;
    double *sigm;
}DATA;

int main_power(DATA* data)
{
    double power,res,kstart,kend,dlc,k,R8=8,R,r,Sigma_R,rho_c,rho_m,rho_mz,H_0=0.1;
    double DZ=1.,Hubble_a,Omega_a;
    int i,n;
    FILE *fd,*fp;
    

//sigma(R)
    rho_c    = 3*(Hubble*Hubble)/(8*PI*G);
    rho_m    = rho_c*Omega;
    
    for(i=0;i<data->n;i++)
    //for(i=0;i<20;i++)
    {
        data->DZ[i]   = GF(1.0/(data->z1[i] + 1), 1.0);
        r = pow(3*data->M[i]/(4*PI)/rho_m,1./3);
        data->sigm[i] = pow(TopHatSigma2(r),0.5);
        printf("DZ,R,M,Sigma_R,rho_m=%lg,%lg,%lg,%lg,%lg\n",data->DZ[i],r/1000,data->M[i],data->sigm[i],rho_m); fflush(stdout);
    }
    printf("end\n");
    return 0;
}

/*    
    char file[512];
    sprintf(file,"data/z%lg-ww%lg-sigma%lg-om%lg-test.dat",redshift,ww,Sigma8,Omega);
    myfopen(fp,file,"w");
    
    for(R = 0.1; R < 200.0; R *= 1.02)
    {
        r = R*1000.;                                                 //unit in kpc/h: *1000
        Sigma_R = pow(TopHatSigma2(r),0.5);
        M = 4*PI*(r*r*r)*rho_m/3;
        fprintf(fp,"%12g%12g%12g%12g%12g\n",DZ,R,M, Sigma_R/DZ,rho_m);
        ///printf("%12g%12g%12g%12g%12g\n",DZ,R,M, Sigma_R/DZ,rho_m);
    }
    fclose(fp);
    printf("%12g%12g\n",R8*1000, pow(TopHatSigma2(R8*1000),0.5));
*/    



