# clean previous data
rm short.dat
rm long.dat
rm etot.dat

N=35

for ss in `seq 7.5 0.05 8.5`
do

LOGFIL=fort.log.$N_$ss

./eggbox_grid_cube.x $N ../../structures/H.xyz $ss 1.0 | tee $LOGFIL

#
# Long part (should be the same for gric cube (a) and ordinary coarse method
#
str=`grep "Ps short a" $LOGFIL`
Ps_short_a=`echo $str | awk '{split($0, a); print a[4]}'`

str=`grep "Ps short t" $LOGFIL`
Ps_short_t=`echo $str | awk '{split($0, a); print a[4]}'`

echo $ss $Ps_short_a $Ps_short_t >> short.dat

#
# Long part (should be the same for gric cube (a) and ordinary coarse method
#
str=`grep "Ps long" $LOGFIL`
Ps_long=`echo $str | awk '{split($0, a); print a[4]}'`

echo $ss $Ps_long >> long.dat

#
# Combined long part and short part
#
str=`grep "E_ps_loc" $LOGFIL`
etot=`echo $str | awk '{split($0, a); print a[3]}'`

str=`grep "Ps total a" $LOGFIL`
etot_a=`echo $str | awk '{split($0, a); print a[4]}'`

echo $ss $etot_a $etot >> etot.dat


done
