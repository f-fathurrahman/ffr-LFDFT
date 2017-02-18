PROGRAM test_fft_fftw3

  IMPLICIT NONE 
  INTEGER :: Nx, Ny, Nz, Npoints
  COMPLEX(8), ALLOCATABLE :: zdata(:)
  INTEGER :: ip

  Nx = 3
  Ny = 5
  Nz = 7
  Npoints = Nx * Ny * Nz

  ALLOCATE( zdata(Npoints) )

  zdata(:) = cmplx(1.d0,1.d0, kind=8)
  zdata(3) = cmplx(2.1,11.7d0, kind=8)

  WRITE(*,'(/,1x,A)') 'Before FFT'
  DO ip = 1, Npoints
    WRITE(*,'(1x,I5,2F18.10)') ip, zdata(ip)
  ENDDO

  CALL fft_fftw3( zdata, Nx, Ny, Nz, .false. )

  WRITE(*,'(/,1x,A)') 'After forward FFT'
  DO ip = 1, Npoints
    WRITE(*,'(1x,I5,2F18.10)') ip, zdata(ip)
  ENDDO

  CALL fft_fftw3( zdata, Nx, Ny, Nz, .true. )

  WRITE(*,'(/,1x,A)') 'After backward FFT'
  DO ip = 1, Npoints
    WRITE(*,'(1x,I5,2F18.10)') ip, zdata(ip)
  ENDDO

  DEALLOCATE( zdata )

  WRITE(*,*) 'Program ended'

END PROGRAM 

