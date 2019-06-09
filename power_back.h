
#define UnitLength_in_cm 3.085678e21   // defines length unit of output (in cm/h)
#define UnitMass_in_g    1.989e43      // defines mass unit of output (in g/cm)
#define UnitVelocity_in_cm_per_s  1e5  // defines velocity unit of output (in cm/sec)
#define PrimordialIndex  0.968
#define HubbleParam      0.71
#define Omega            0.275
#define OmegaLambda      0.725         // Cosmological constant (at z=0)
#define OmegaBaryon      0.045
#define PI               3.14159265358979323846
#define  GRAVITY         6.672e-8
#define  HUBBLE          3.2407789e-18   /* in h/sec */
#define Sigma8           0.816//0.85//0.63297
#define ww               -1.
double G, Hubble,UnitTime_in_s;


extern double PowerSpec_DM_2ndSpecies(double k);
extern double PowerSpec_EH(double k);
extern double tk_eh(double k);
extern double TopHatSigma2(double R);
extern double sigma2_int(double k);
extern void set_units(void);
void initialize_powerspectrum(void);
