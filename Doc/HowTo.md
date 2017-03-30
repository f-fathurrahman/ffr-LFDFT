

## Initializing LF grids

```fortran
! Define the box
AA = (/ 0.d0, 0.d0, 0.d0 /)
BB = (/ 6.d0, 6.d0, 6.d0 /)
! define number of sampling points
NN = (/ 25, 25, 25 /)
! This will initialize global variables at `m_LF3d`
CALL init_LF3d_p( NN, AA, BB )
```

## Using 3D FFT from FFTW3

```fortran
CALL fft_fftw3( zdata, Nx, Ny, Nz, .FALSE. )  ! forward transform

call fft_fftw3( zdata, Nx, Ny, Nz, .TRUE. ) ! inverse transform
```

## Setting up calculation

Allocate memory for Hamiltonian (potentials and electron density)

```fortran
CALL allocate_hamiltonian()
```

Set up local potential (harmonic potential, for example):

```fortran
CALL init_V_pc_loc_harmonic( omega, center )
```

Initialize electronic states and their occupations and also
arrays for Kohn-Sham eigenvectors (orbitals) and eigenvalues:

```fortran
! Initialize electronic states variables
Nstates = 4
ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )
ALLOCATE( Focc(Nstates) )
Focc(:) = 2.d0
```

Random initialization of KS eigenvectors:
```fortran
DO ist = 1, Nstates
  DO ip = 1, Npoints
    CALL random_number( evecs(ip,ist) )
  ENDDO
ENDDO
CALL orthonormalize( Nstates, evecs )
```

Solve Kohn-Sham equations via direct minimization (conjugate gradient
algorithm)

```
CALL KS_solve_Emin_pcg( 3.d-5, 100, .FALSE. )
```

Free memory

```
CALL dealloc_hamiltonian()
CALL dealloc_LF3d()
```
