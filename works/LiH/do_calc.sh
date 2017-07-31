for nn in 11 13 15 17 19 21 23 27 33 35 39 41 45 49 55
do

cat << EOF > INPUT
&CONTROL
  pseudo_dir = '../../HGH'
  etot_conv_thr = 1.0d-6
/

&SYSTEM
  ibrav = 8
  nat = 2
  ntyp = 2
  A = 8.46683536902
  B = 8.46683536902
  C = 8.46683536902
  nr1 = $nn
  nr2 = $nn
  nr3 = $nn
/

&ELECTRONS
  KS_Solve = 'Emin_pcg'
  cg_beta = 'PR'
  electron_maxstep = 150
  mixing_beta = 0.5
  diagonalization = 'LOBPCG'
/

ATOMIC_SPECIES
Li   3.0  Li_sc.hgh
H    1.0  H.hgh

ATOMIC_POSITIONS angstrom
Li      4.23341768       4.23341768       5.04089768
H       4.23341768       4.23341768       3.42593768

EOF

LOGFILE="fort.log.$nn"
../../ffr_LFDFT_pgi.x INPUT | tee $LOGFILE

str=`grep "! Electronic" $LOGFILE`
etot=`echo $str | awk '{split($0, a); print a[4]}'`

str=`grep spacing $LOGFILE`
spacing=`echo $str | awk '{split($0, a); print a[4]}'`

echo $spacing $etot >> RESULT.dat

done
