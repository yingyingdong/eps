FILE *logfile;// logfile=stdout; if needed.
#define myfopen(filepointer,filename,filemode) if(!((filepointer)=fopen(filename,filemode))){ fprintf(logfile,"Error opening file '%s'\n",filename);    fflush(logfile); exit(1);   }
double GF(double astart, double aend);
double EZ_dfy(double a);
double hubble_dfy(double a);
double Omega_dfy(double a);
