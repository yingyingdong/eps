#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Romberg.h"
#include "power.h"
#include "growth.h"

double M2R(double M,double rho_m)
{
    if(tophat==1) return pow(3*M/(4*PI)/rho_m,1./3);
    else          return pow(M/(6*PI*PI)/rho_m,1./3);
}

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

int main_power(int n, double *zz, double *M, double *DZ, double *sigm)
{
    double R8=8,R,r,Sigma_R,rho_c,rho_m,rho_mz,H_0=0.1;
    double Hubble_a,Omega_a;
    int i;
    FILE *fd,*fp;
    

//sigma(R)
    rho_c    = 3*(Hubble*Hubble)/(8*PI*G);
    rho_m    = rho_c*Omega;
    
    for(i=0;i<n;i++)
    //for(i=0;i<20;i++)
    {
        DZ[i]   = GF(1.0/(zz[i] + 1), 1.0);
        r       = M2R(M[i],rho_m);
        sigm[i] = pow(TopHatSigma2(r),0.5);
        printf("DZ,R,M,Sigma_R,rho_m=%lg,%lg,%lg,%lg,%lg\n",DZ[i],r/1000,M[i],sigm[i],rho_m); fflush(stdout);
    }
    printf("end\n");
    return 0;
}

int dc_redshift(double *dc,double *zz,int nz)
{
    int i;
    double DZ;

    for(i=0;i<nz;i++)
    {
        DZ    = GF(1.0/(zz[i] + 1), 1.0);
        dc[i] = 1.686*DZ;
    }

    return 0;
}

int sig_M1_M2_zz(double M, double *dM, double *S2,int n, double zz, double *temp)
{
    int i;
    double DZ,rho_c,rho_m,r1,r2,S1;
    //S2    = v=(double *)malloc((size_t) ((n)*sizeof(double))); 
    
    temp[1] = GF(1.0/(zz + 1), 1.0);
    rho_c   = 3*(Hubble*Hubble)/(8*PI*G);
    rho_m   = rho_c*Omega;
    
    r1      = M2R(M,rho_m);
    temp[0] = TopHatSigma2(r1); 

    for(i=0;i<n;i++)
    {
        if(tophat==1) r2 = M2R(M+dM[i],rho_m);
        S2[i] = TopHatSigma2(r2); //sigm^2
    }
    return 0;
}

int sig_M1_M2_z1_z2(double M1, double *M2, double *S2,int n, double z1, double z2,double *temp)
{
    int i;
    double DZ,rho_c,rho_m,r1,r2,S1;
    
    temp[1] = GF(1.0/(z1 + 1), 1.0);
    temp[2] = GF(1.0/(z2 + 1), 1.0);
    rho_c   = 3*(Hubble*Hubble)/(8*PI*G);
    rho_m   = rho_c*Omega;

    r1      = M2R(M1,rho_m);
    temp[0] = TopHatSigma2(r1);

    for(i=0;i<n;i++)
    {
        if(tophat==1) r2 = M2R(M2[i],rho_m);
        S2[i] = TopHatSigma2(r2); //sigm^2
    }
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



