#!/bin/bash
for i in `seq 1 20`
do

str=`grep "Total      =" LOG_*`
E_ps_loc_ss=`echo $str | awk '{split($0, a); print a[3]}'`

echo

done
