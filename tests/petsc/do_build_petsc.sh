#!/bin/bash

INC="-I/home/efefer/mysoftwares/petsc-3.7.5/include/ -I../../"
LIB="../../libmain.a libLFDFT_petsc.a -lblas -llapack -lfftw3 \
-L/home/efefer/mysoftwares/petsc-3.7.5/lib -lpetsc"

make

bas=`basename $1 .F`
# remove the previous executable
rm -vf $bas.x
mpifort -Wall -O3 -ffree-form $INC $1 $LIB -o $bas.x

# for
#mpifort -free $INC $1 $LIB -o $bas.x
