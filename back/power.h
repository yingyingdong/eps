
//#define UnitLength_in_cm 3.085678e24   // defines length unit of output (in cm/h)
#define UnitLength_in_cm 3.085678e21   // defines length unit of output (in cm/h)
#define UnitMass_in_g    1.989e43      // defines mass unit of output (in g/cm)
#define UnitVelocity_in_cm_per_s  1e5  // defines velocity unit of output (in cm/sec)
#define PrimordialIndex  1.//0.9645//1.
#define HubbleParam      0.71//0.6727//0.71//0.679//0.71
#define Omega            0.268//0.3156//0.3//0.268//0.2568 //0.268//0.2568//0.268
#define OmegaLambda      0.732//0.6844//0.7//0.732//0.7432 //0.732//0.7432//0.732         // Cosmological constant (at z=0)
#define OmegaBaryon      0.045//0.049168522713946//0.045
#define PI               3.14159265358979323846
#define  GRAVITY         6.672e-8
#define  HUBBLE          3.2407789e-18   /* in h/sec */
#define Sigma8           0.85//0.831 //1.1414//0.85//0.8429//0.85//0.816//0.85//0.63297
#define ww               -1.
double G, Hubble,UnitTime_in_s;


extern double PowerSpec_DM_2ndSpecies(double k);
extern double PowerSpec_EH(double k);
extern double tk_eh(double k);
extern double TopHatSigma2(double R);
extern double sigma2_int(double k);
extern void set_units(void);
void initialize_powerspectrum(void);
