#/bin/bash
for nn in 11 13 15 17 19 21 23 27 33 35 39 41 45 49 55
do
  ./do_Emin_pcg_gaussian_G.x $nn $1 $2 | tee fort.log.$nn
done

grep "! Elec" fort.log.* > RESULT.dat

