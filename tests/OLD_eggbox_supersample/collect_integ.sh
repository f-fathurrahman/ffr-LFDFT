# force to use dot for decimal separator
export LC_NUMERIC=en_US.UTF-8
element=$1
center=`seq 0.0 0.1 16.0`
Npoints="55"


for c in $center
do
  echo "center = $c"
  echo "1" > ATOM.xyz
  echo "" >> ATOM.xyz
  echo "$element 0.0 $c 0.0" >> ATOM.xyz
  for N in $Npoints
  do
    ./test_eggbox.x $N
  done
done

