# These files are from the CERNLIB package. To compile a library file (using gfortran), first make the object files:
gfortran -c -fPIC wclbes.F clogam64.F r8dp.F r7dp.F r6dp.F r5dp.F r4dp.F r3dp.F r2dp.F r1dp.F cdigam64.F mtlset.F lenocc.F abend.F -D CERNLIB_DOUBLE
# Next make the .a archive file from all the relevant object files
ar r libwclbes.a wclbes.o clogam64.o r8dp.o r7dp.o r6dp.o r5dp.o r4dp.o r3dp.o r2dp.o r1dp.o cdigam64.o mtlset.o mtlprt.o lenocc.o abend.o
# Finally, link the ckin to this archive e.g.:
gcc ckin.c -Wall -o ckin -lm -O2 -DHAVE_WCLBES -L/home/path/to/wclbes -lwclbes -lgfortran

#To compile a library file (using f77/f90), first make the object files:
f77 -c -fPIC wclbes.F clogam64.F r8dp.F r7dp.F r6dp.F r5dp.F r4dp.F r3dp.F r2dp.F r1dp.F cdigam64.F mtlset.F mtlprt.F lenocc.F abend.F -D CERNLIB_DOUBLE

# Next make the .a archive file from all the relevant object files
ar r libwclbes.a wclbes.o clogam64.o r8dp.o r7dp.o r6dp.o r5dp.o r4dp.o r3dp.o r2dp.o r1dp.o cdigam64.o mtlset.o mtlprt.o lenocc.o abend.o

# Finally, link the program to these. The program will need to be compiled first and then ALWAYS linked using f77/f90 e.g.:
gcc -c ckin.c -o ckin.o -DHAVE_WCLBES
f90 -o ckin ckin.o -L/export/home/user/software/wclbes -lwclbes -lm -lfsu
