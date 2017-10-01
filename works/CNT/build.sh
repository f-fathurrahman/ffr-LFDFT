#!/bin/bash

INC="-I../../"
LIB="../../libmain.a -lblas -llapack -lfftw3"

bas="CNT_main"

# remove the previous executable
rm -vf $bas.x

$2 $INC $1 $LIB -o $bas.x

echo "Test executable: $bas.x"

