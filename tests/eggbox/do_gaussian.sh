NPOINTS=`seq 11 2 45`
for n in $NPOINTS
do
  ./test_gaussian_G.x $n 2.5 & python plot_cols.py $n
done

