# ffr-LFDFT

This is a work in progress.

## Introduction

An experimental package to solve electronic structure based on density functional theory
using Lagrange basis functions.

It can do total energy calculations via SCF and direct minimization.

Tested compilers:
- [`gfortran`](https://gcc.gnu.org/fortran/)
- [`g95`](www.g95.org)
- [`ifort`](https://software.intel.com/en-us/fortran-compilers)
- [`pgf90`](https://www.pgroup.com/products/community.htm)
- [`sunf95`](http://www.oracle.com/technetwork/server-storage/developerstudio/downloads/index.html)

Dependencies
- BLAS and LAPACK
- FFTW3
- LibXC

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


