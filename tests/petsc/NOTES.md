## PETSC notes

I used petsc-3.7.5 and installed two versions of it: `debug` and `opt`.
`debug` is the default version that I used for development.

Here is a dumb script to compile and link the code
```bash
INC="-I/home/efefer/mysoftwares/petsc-3.7.5/include/ -I../../"
LIB="../../libmain.a -lblas -llapack -lfftw3 \
-L/home/efefer/mysoftwares/petsc-3.7.5/lib -lpetsc"

bas=`basename $1 .F`
# remove the previous executable
rm -vf $bas.x
mpifort -Wall -O3 -ffree-form $INC $1 $LIB -o $bas.x
```

I originally thought that I need to use different `INC` and `LIB`
if I want to test the `opt` version. It seems that there is I don't need
to do that.
I tried to do that before, however, because I haven't added the
`LD_LIBRARY_PATH` for the `opt`, I always get the `debug` version,
even if I passed the path to `opt` version of the library in the
compile-link process.
Simply appending `opt` library path to `LD_LIBRARY_PATH` will enable
the of `opt` version of the library.
