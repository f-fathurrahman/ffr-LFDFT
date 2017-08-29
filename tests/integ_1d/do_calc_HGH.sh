for N in `seq 3 2 55`; do
  for scal in `seq 0.0 0.1 0.9`; do
    ./HGH_integ_1d.x $N $scal ../../HGH/H.hgh
  done
  echo "#"
done

