for N in `seq 15 4 75`; do
  rm res_$N.dat
  for scal in `seq 0.0 0.025 1.0`; do
    ./HGH_integ_1d.x $N $scal ../../HGH/H.hgh >> res_$N.dat
  done
done

