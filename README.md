# ffr-LFDFT

This is a work in progress.

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

Main executable is named `ffr_LFDFT_<compiler_name>.x`, for example `ffr_LFDFT_gfortran.x`

The input files are very similar to Quantum Espresso's `pw.x` input.

To run the main program

```
ffr_LFDFT_gfortran.x INPUT > LOG
```

One post-processing program, with very limited capability, is also provided.
Currently it only can produce 3D orbital in XSF (Xcrysden) format.


## Examples

Examples can be found in `work` directory.


