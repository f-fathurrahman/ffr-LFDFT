#!/bin/bash

INC="-I../../src/"
LIB="../../src/libmain.a ../../src/libpoisson_ISF.a ../../src/libsparskit.a -lblas -llapack -lfftw3 \
-L/home/efefer/mysoftwares/libxc-3.0.0/lib -lxcf90 -lxc"

bas=`basename $1 .f90`

# remove the previous executable
rm -vf $bas.x

$2 $INC $1 $LIB -o $bas.x

echo "Test executable: $bas.x"

