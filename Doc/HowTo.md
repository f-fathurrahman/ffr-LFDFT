

## Initializing LF grids

Three types of `LF3d`:

- periodic LF, initialized using `init_LF3d_p()`
- box/cluster LF, initialized using `init_LF3d_c()`
- sinc LF, initialized using `init_LF3d_sinc()`

Example code for constructing periodic LF:

```fortran
! Define the box
AA = (/ 0.d0, 0.d0, 0.d0 /) ! left boundaries
BB = (/ 6.d0, 6.d0, 6.d0 /) ! right boundaries
NN = (/ 25, 25, 25 /) ! define number of sampling points
CALL init_LF3d_p( NN, AA, BB ) ! This will initialize global variables at `m_LF3d`
```

Similiarly for box/cluster LF:
```fortran
CALL init_LF3d_c( NN, AA, BB )
```

Meanwhile, for sinc LF:
```fortran
hh(1:3) = BB(:) - AA(:) / (NN(:) - 1)
CALL init_LF3d_sinc( NN, hh )
```

## Getting arguments

```fortran
CHARACTER(56) :: args_chars
INTEGER :: iargc
INTEGER :: N
CHARACTER(56) :: fname

!
IF( iargc() /= 1 ) THEN
  WRITE(*,*) 'Need exactly two arguments'
  STOP
ENDIF

! Get first argument
CALL getarg( 1, args_chars )
fname = args_chars

! Get second argument
CALL getarg( 2, args_chars)
! Convert string to integer
READ( args_chars, * ) N
```

## Using 3D FFT from FFTW3

This is important for Poisson solver for periodic BC:

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

```fortran
CALL KS_solve_Emin_pcg( 3.d-5, 100, .FALSE. )
```

Free memory

```fortran
CALL dealloc_hamiltonian()
CALL dealloc_LF3d()
```
