#/bin/bash
rm -fv PSolver.o
g95 -fno-second-underscore -I../../ MAIN.f90 init_density.f90 ../../libmain.a *.o ABINIT-common/*.o
