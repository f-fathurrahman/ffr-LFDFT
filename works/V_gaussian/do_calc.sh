for nn in 11 13 15 17 19 21 23 27 33 35
do
  ./do_Emin_pcg_gaussian_G.x $nn > fort.log.$nn
done

grep "Total    =" fort.log.* > RESULT.dat

