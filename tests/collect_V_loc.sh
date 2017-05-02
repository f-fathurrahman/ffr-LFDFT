element=$1
center="0.0 8.0  7.0"
Npoints="55 65 75 85"


for c in $center
do
  echo "center = $c"
  echo "1" > STRUCT.xyz
  echo "" >> STRUCT.xyz
  echo "$element 0.0 $c 0.0" >> STRUCT.xyz
  for N in $Npoints
  do
    ./do_print_V_ps_loc.x $N STRUCT.xyz
    mv fort.$N fort.$N.$element.$c.shifted
  done
done
