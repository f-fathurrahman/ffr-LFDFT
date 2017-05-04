!! PURPOSE:
!!
!!   This subroutine calculates second derivative
!!   matrix of sinc Lagrange basis functions in 1D.
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! NOTE:
!!
!!   This subroutine uses full storage of the matrix,
!!   not yet exploiting the symmetry of derivative matrix.
!!
!! IMPORTANT:
!!
!!   The matrices D2jl are assumed to have been allocated somewhere else.
!!

SUBROUTINE init_deriv_matrix_sinc( N, h, grid, D2jl )
  USE m_constants, ONLY : PI
  IMPLICIT NONE

  !! Number of points
  INTEGER :: N

  !! Scaling factor
  REAL(8) :: h

  !! Grid points
  REAL(8) :: grid(N)

  !! Second derivative matrix
  REAL(8) :: D2jl(N,N)

  INTEGER :: i, j
  
  ! Diagonal part
  DO i = 1, N
    D2jl(i,i) = -pi**2/3.d0/h**2
  ENDDO
  ! Off-diagonal
  DO j=1,N
    DO i=j+1,N
      D2jl(i,j) = -2.d0*(-1)**(i-j)/( grid(i)-grid(j) )**2
      D2jl(j,i) = D2jl(i,j)
    ENDDO
  ENDDO

END SUBROUTINE 
