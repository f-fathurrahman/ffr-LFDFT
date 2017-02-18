for N in `seq 17 2 41`
do
  echo
  echo 'Testing pcg:' $N 
  ./a.out $N pcg
done
