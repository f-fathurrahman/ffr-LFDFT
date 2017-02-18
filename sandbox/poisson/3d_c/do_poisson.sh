for N in `seq 17 2 41`
do
  echo
  echo 'Testing cg:' $N 
  ./p2.x $N cg
done
