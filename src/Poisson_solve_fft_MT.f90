!! PURPOSE
!!
!!   This subroutine solves Poisson equation using Fast Fourier
!!   transform.
!!
!! AUTHOR
!!
!!   Fadjar Fathurrahman
!!
!! NOTES
!!
!!   The input `rho` will be multiplied by -4*pi.
!!   The output is given in `phi`.

SUBROUTINE Poisson_solve_fft_MT( rho, phi )
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     LL => LF3d_LL, &
                     NN => LF3d_NN, &
                     G2 => LF3d_G2
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: rho(Npoints)
  REAL(8) :: phi(Npoints)
  ! Local
  COMPLEX(8), ALLOCATABLE :: tmp_fft(:)
  INTEGER :: ip, Nx, Ny, Nz
  REAL(8) :: R, RG

  ! shortcuts
  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  R = minval( LL )*0.5d0
  !WRITE(*,*) 'Poisson_solve_fft_MT: R = ', R

  ALLOCATE( tmp_fft(Npoints) )
  DO ip = 1, Npoints
    tmp_fft(ip) = cmplx( rho(ip), 0.d0, kind=8 )
  ENDDO

  ! forward FFT
  CALL fft_fftw3( tmp_fft, Nx, Ny, Nz, .false. )  ! now `tmp_fft = rho(G)`
  
  tmp_fft(1) = 2.d0*PI * R**2 * tmp_fft(1)  ! V_H(G=0)

  DO ip = 2, Npoints
    RG = R*sqrt(G2(ip))
    tmp_fft(ip) = 4.d0*PI*tmp_fft(ip) / G2(ip) * (1.d0-cos(RG))
  ENDDO

  ! Inverse FFT
  CALL fft_fftw3( tmp_fft, Nx, Ny, Nz, .true. )

  DO ip = 1, Npoints
    phi(ip) = real( tmp_fft(ip), kind=8 )
  ENDDO

  DEALLOCATE( tmp_fft )

END SUBROUTINE

