for N in `seq 15 4 75`; do
  rm res_$N.dat
  for scal in `seq 0.0 0.05 0.9`; do
    ./HGH_integ_1d.x $N $scal ../../HGH/H.hgh >> res_$N.dat
  done
done

