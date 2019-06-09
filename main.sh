#gcc -c Romberg.c power.c growth.c main_power.c -lm -I.    #for power.c
#ar -crv ps.a *.o #static library

\rm *.o *.so
gcc -c -fPIC Romberg.c power.c growth.c main_power.c -lm -I.
gcc -shared -o pshmf.so *.o

