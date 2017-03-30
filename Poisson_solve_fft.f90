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
!!   This subroutine should be called for periodic system only.

SUBROUTINE Poisson_solve_fft( rho, phi )
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     NN => LF3d_NN, &
                     G2 => LF3d_G2
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: rho(Npoints)
  REAL(8) :: phi(Npoints)
  ! Local
  COMPLEX(8), ALLOCATABLE :: tmp_fft(:)
  INTEGER :: ip, Nx, Ny, Nz

  ! shortcuts
  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  ALLOCATE( tmp_fft(Npoints) )
  DO ip = 1, Npoints
    tmp_fft(ip) = cmplx( rho(ip), 0.d0, kind=8 )
  ENDDO

  ! forward FFT
  CALL fft_fftw3( tmp_fft, Nx, Ny, Nz, .false. )  ! now `tmp_fft = rho(G)`
  tmp_fft(1) = (0.d0,0.d0)
  DO ip = 2, Npoints
    tmp_fft(ip) = 4.d0*PI*tmp_fft(ip) / G2(ip)
  ENDDO
  ! now `tmp_fft` = phi(G)

  ! Inverse FFT
  CALL fft_fftw3( tmp_fft, Nx, Ny, Nz, .true. )

  DO ip = 1, Npoints
    phi(ip) = real( tmp_fft(ip), kind=8 )
  ENDDO

  DEALLOCATE( tmp_fft )


END SUBROUTINE

