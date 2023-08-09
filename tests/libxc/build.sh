#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo
  echo "ERROR"
  echo "Need two parameters: main file and compiler name (and or options)"
  echo "Example: ./build.sh myfile.f90 \"pgf90 -O2\""
  echo
  exit 1
fi

LIBXC_HOME=/home/efefer/mysoftwares/libxc-5.1.5
LIB="-L$LIBXC_HOME/lib -lxcf90 -lxc"

bas=`basename $1 .f90`

# remove the previous executable
rm -vf $bas.x

$2 libxc_funcs.f90 libxc.f90 $1 $LIB -o $bas.x

#gfortran -Wall -O3 -ffree-form $INC $1 $LIB -o $bas.x
echo "Test executable: $bas.x"

