N=45
Rcut="1.5"

for ss in `seq 7.5 0.05 8.5`
#for ss in "8.0"
do

LOGFIL=fort.log.$N_$ss
LOGFIL_STD=fort.log.$N_${ss}_STD

./test1.x $N ../../HGH/H.hgh 8.0 $ss $Rcut > $LOGFIL

str=`grep "E_ps_loc" $LOGFIL`
E_ps_loc_ss=`echo $str | awk '{split($0, a); print a[3]}'`

# Standard algorithm
./standard.x $N ../../HGH/H.hgh 8.0 $ss $Rcut > $LOGFIL_STD
str=`grep "E_ps_loc" $LOGFIL_STD`
E_ps_loc_std=`echo $str | awk '{split($0, a); print a[3]}'`

echo $ss $E_ps_loc_ss $E_ps_loc_std >> Ene_N_${N}_rcut_${Rcut}.dat


done
