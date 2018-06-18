#/bin/bash
#for nn in 10 12 14 16 18 20 22 26 32 34 38 40 44 48 54
#do
#  ./do_Emin_pcg_gaussian_sinc.x $nn $1 $2 | tee fort.log.$nn
#done

grep "! Elec" fort.log.* > RESULT.dat

