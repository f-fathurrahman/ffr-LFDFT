

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

CALL fft_fftw3( zdata, Nx, Ny, Nz, .TRUE. ) ! inverse transform
```

## Interpolation using `bspline.f90` for periodic function

We have data `ctmp(Nx*Ny*Nz)`, defined on FFT grid (only one end-points on FFT
box).

Grid points for periodic LF are different from the usual FFT grid.
We need to calculate / shift the grid points

```fortran
shiftx = 0.5d0*( grid_x(2) - grid_x(1) )
DO i = 1,Nx
  x(i) = grid_x(i) - shiftx
ENDDO
x(Nx+1) = x(Nx) + 2*shiftx

shifty = 0.5d0*( grid_y(2) - grid_y(1) )
DO j = 1,Ny
  y(j) = grid_y(j) - shifty
ENDDO
y(Ny+1) = y(Ny) + 2*shifty

shiftz = 0.5d0*( grid_z(2) - grid_z(1) )
DO k = 1,Nz
  z(k) = grid_z(k) - shiftz
ENDDO
z(Nz+1) = z(Nz) + 2*shiftz
```

We also need to copy the linearized version of our 3D data to appropriate
3D array (FFT grid points + padded boundary points):

```fortran
DO kk = 1,Nz+1
DO jj = 1,Ny+1
DO ii = 1,Nx+1
  i = ii
  j = jj
  k = kk
  IF( kk == Nz+1 ) k = 1
  IF( jj == Ny+1 ) j = 1
  IF( ii == Nx+1 ) i = 1
  ip = xyz2lin(i,j,k)
  !
  data3d(ii,jj,kk) = real( ctmp(ip), kind=8 )
ENDDO
ENDDO
ENDDO
```

These variables are needed for initialization of splines:

```fortran
iknot = 0  ! determine the knots automatically
kx = 4
ky = 4
kz = 4
ALLOCATE( tx(Nx+1+kz), ty(Ny+1+ky), tz(Nz+1+kz) )
ALLOCATE( x(Nx+1), y(Ny+1), z(Nz+1) )
```

Additional variables (**TODO**):

```fortran
iloy = 1
iloz = 1
inbvx = 1
inbvy = 1
inbvz = 1
idx = 0
idy = 0
idz = 0
```


We can then initialize the splines:

```fortran
CALL db3ink( x, Nx+1, y, Ny+1, z, Nz+1, data3d, kx,ky,kz, iknot, tx,ty,tz, data3d, iflag )
```

```fortran
CALL db3val( dx, dy, dz, idx,idy,idz, tx,ty,tz, Nx+1,Ny+1,Nz+1,kx,ky,kz, data3d, &
       val, iflag, inbvx, inbvy, inbvz, iloy, iloz )
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
CALL KS_solve_Emin_pcg( 3.d-5, .FALSE. )
```

Free memory

```fortran
CALL dealloc_hamiltonian()
CALL dealloc_LF3d()
```
