element=$1
center="0.0 0.5 1.0 1.5 2.0"
Npoints="35"


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

