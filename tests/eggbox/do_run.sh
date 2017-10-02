rm -f RESULT.dat
for xx in `seq 8.0 0.1 10.0`
do
cat << EOF > ATOM.xyz
1

Be_sc  $xx 8.0 8.0
EOF

res=`./test_grid_atom_cube.x 45`
echo $xx $res

done
