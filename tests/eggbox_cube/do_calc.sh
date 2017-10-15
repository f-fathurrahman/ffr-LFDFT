N=45
for ss in `seq 7.5 0.1 8.5`
do

LOGFIL=fort.log.$N_$ss
#./eggbox_grid_cube.x $N ../../structures/H.xyz $ss | tee $LOGFIL

str=`grep "Ps short a" $LOGFIL`
Ps_short_a=`echo $str | awk '{split($0, a); print a[4]}'`

str=`grep "Ps short t" $LOGFIL`
Ps_short_t=`echo $str | awk '{split($0, a); print a[4]}'`

str=`grep "Ps long" $LOGFIL`
Ps_long=`echo $str | awk '{split($0, a); print a[4]}'`

#echo $ss $Ps_short_a $Ps_short_t >> short.dat
#echo $ss $etot >> RESULT.dat

echo $ss $Ps_long >> long.dat

done
