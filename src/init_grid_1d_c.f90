!! PURPOSE:
!!
!!   This subroutine initializes cluster Lagrange basis functions
!!   in 1D, defined on [A,B] with N points.
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! IMPORTANT:
!!
!!   The array `grid` is assumed to be already allocated properly elsewhere

SUBROUTINE init_grid_1d_c( N, A, B, grid )

  IMPLICIT NONE
  ! Arguments
  INTEGER :: N
  REAL(8) :: A, B
  REAL(8) :: grid(N)
  ! Local
  INTEGER :: ii

  ! No need to check if B > A. This has been done at init_LF3d_c

  ! Initialization of grid points
  DO ii=1,N
    grid(ii) = A + ii*(B-A)/(N+1)
  ENDDO

END SUBROUTINE

