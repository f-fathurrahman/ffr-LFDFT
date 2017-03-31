!! PURPOSE:
!!
!!   This subroutine calculates first derivative matrix D1jl and second derivative
!!   matrix of Lagrange basis functions in 1D.
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! NOTE:
!!
!!   This subroutine uses full storage of the matrix,
!!   not yet exploiting the symmetry of D1jl and D2jl
!!
!! IMPORTANT:
!!
!!   The matrices D1jl and D2jl are assumed to have been allocated somewhere else.
!!

SUBROUTINE init_deriv_matrix_c( N, L, D2jl )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  
  !! Number of points
  INTEGER :: N
  
  !! Length of domain
  REAL(8) :: L

  !! Second derivative matrix
  REAL(8) :: D2jl(N,N)

  ! Local
  REAL(8) :: t1, t2, pre
  INTEGER :: nnm, nnp
  INTEGER :: ll, jj

  ! Diagonal
  pre = -PI**2/(2*L**2)
  DO ll = 1, N
    t1 = ( 2*(N+1)**2 + 1 )/3.d0
    t2 = sin( PI*ll/(N+1) )**2
    D2jl(ll,ll) = pre*(t1 - 1.0/t2)
  ENDDO

  ! Off-diagonal
  DO ll = 1, N
    DO jj = ll+1, N
      nnm = ll - jj
      nnp = ll + jj
      pre = -PI**2/(2*L**2)*(-1)**nnm
      t1 = sin( PI*nnm/2.d0/(N+1) )**2
      t2 = sin( PI*nnp/2.d0/(N+1) )**2
      !
      D2jl(ll,jj) = pre*(1.d0/t1 - 1.d0/t2)
      D2jl(jj,ll) = pre*(1.d0/t1 - 1.d0/t2)
    ENDDO
  ENDDO

END SUBROUTINE 
