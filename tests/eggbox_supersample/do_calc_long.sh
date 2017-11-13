N=55
Rcut="1.0"  # not really used in the calculation

for ss in `seq 7.5 0.05 8.5`
do

LOGFIL_LONG=fort.log.${N}_LONG

# long part only
./test_V_long_only.x $N ../../HGH/H.hgh 8.0 $ss $Rcut > $LOGFIL_LONG

str=`grep "E_ps_loc" $LOGFIL_LONG`
E_ps_loc_long=`echo $str | awk '{split($0, a); print a[3]}'`

echo $ss $E_ps_loc_long >> Ene_N_${N}_long.dat

done
