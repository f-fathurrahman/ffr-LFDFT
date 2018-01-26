for N in `seq 3 2 55`; do
  for scal in `seq 0.0 0.1 0.9`; do
    ./integ_1d.x $N $scal 1.0 5.5
  done
  echo "#"
done

