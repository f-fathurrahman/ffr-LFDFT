N=45
Rcut="1.0"  # dummy, not really used in the calculation
ELEM="C"

for ss in `seq 7.5 0.05 8.5`
do

LOGFIL_STD=fort.log.${N}_${ss}_STD

# Standard algorithm
./standard.x $N ../../HGH/${ELEM}.hgh 8.0 $ss $Rcut > $LOGFIL_STD
str=`grep "E_ps_loc" $LOGFIL_STD`
E_ps_loc_std=`echo $str | awk '{split($0, a); print a[3]}'`

echo $ss $E_ps_loc_std >> Ene_N_${N}_std.dat

done

