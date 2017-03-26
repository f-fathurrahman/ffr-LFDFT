#!/bin/bash

INC="-I/home/efefer/mysoftwares/petsc-3.7.5/include/ -I../../"
LIB="../../libmain.a -lblas -llapack -lfftw3 \
-L/home/efefer/mysoftwares/petsc-3.7.5/lib -lpetsc"

bas=`basename $1 .F`

mpifort -Wall -O3 -ffree-form $INC $1 $LIB -o $bas.x

# for
#mpifort -free $INC $1 $LIB -o $bas.x

