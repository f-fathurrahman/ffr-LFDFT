#!/bin/bash

INC="-I../"
LIB_SPARSKIT=/home/efefer/mysoftwares/lib/libskit.a
LIB="../libmain.a -lblas -llapack -lfftw3 $LIB_SPARSKIT"

cd ../
make Makefile
cd tests/

bas=`basename $1 .f90`

# remove the previous executable
rm -vf $bas.x

gfortran -Wall -O3 -ffree-form $INC $1 $LIB -o $bas.x
echo "Test executable: $bas.x"

# for
#mpifort -free $INC $1 $LIB -o $bas.x
