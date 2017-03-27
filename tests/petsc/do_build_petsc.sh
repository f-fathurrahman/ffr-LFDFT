#!/bin/bash

INC="-I/home/efefer/mysoftwares/petsc-3.7.5/include/ -I../../"
LIB="../../libmain.a -lblas -llapack -lfftw3 \
-L/home/efefer/mysoftwares/petsc-3.7.5/lib -lpetsc"

INC_OPT="-I/home/efefer/mysoftwares/petsc-3.7.5_OPT/include/ -I../../"
LIB_OPT="../../libmain.a -lblas -llapack -lfftw3 \
-L/home/efefer/mysoftwares/petsc-3.7.5_OPT/lib -lpetsc"

bas=`basename $1 .F`

# remove the previous executable
rm -vf $bas.x
rm -vf "$bas"_opt.x

mpifort -Wall -O3 -ffree-form $INC $1 $LIB -o $bas.x

mpifort -Wall -O3 -ffree-form $INC_OPT $1 $LIB_OPT -o "$bas"_opt.x

# for
#mpifort -free $INC $1 $LIB -o $bas.x

