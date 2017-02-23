

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


