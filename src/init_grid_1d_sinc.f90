!! PURPOSE:
!!
!!   This subroutine initializes sinc Lagrange basis functions
!!   in 1D, defined with scaling factor h and N points.
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! IMPORTANT:
!!
!!   The array `grid` is assumed to be already allocated properly elsewhere

SUBROUTINE init_grid_1d_sinc( N, h, grid )

  IMPLICIT NONE
  ! Arguments
  INTEGER :: N
  REAL(8) :: h
  REAL(8) :: grid(N)
  ! Local
  INTEGER :: ii
  REAL(8) :: A

  A = -(N-1)/2.d0*h
  ! Initialization of grid points
  DO ii=1,N
    grid(ii) = A + (ii-1)*h
  ENDDO

END SUBROUTINE

