#!/bin/bash

INC="-I../"
LIB="../libmain.a -lblas -llapack -lfftw3"

bas=`basename $1 .f90`

# remove the previous executable
rm -vf $bas.x

$2 $INC $1 $LIB -o $bas.x

#gfortran -Wall -O3 -ffree-form $INC $1 $LIB -o $bas.x
echo "Test executable: $bas.x"

# for
#mpifort -free $INC $1 $LIB -o $bas.x
