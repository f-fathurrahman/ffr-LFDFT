# ffr-LFDFT

## Introduction

An experimental package to solve electronic structure based on density functional theory
using Lagrange basis functions.

Tested compilers:
- gfortran
- g95
- ifort
- pgi
- sunf95

Dependencies
- BLAS and LAPACK
- FFTW3

Contains selected `SPARSKIT` files.

Modified `bspline-fortran`.

## Building

Modify file `make.inc` to suit your needs.

```
cd src
make           # build the library
make main      # for main executable
make postproc  # for post processing
```

## Examples

Examples can be found in `work` directory.

