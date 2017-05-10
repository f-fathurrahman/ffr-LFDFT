! Evaluate the value of L(x) for arbitrary x according to the analytic form
! of the basis function.
!------------------------------------------------------------------------------
FUNCTION eval_LF1d_p(N, L, grid, ibf, x) RESULT(ff)
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  ! arguments
  INTEGER :: N
  REAL(8) :: L
  REAL(8) :: grid(N)
  INTEGER :: ibf
  REAL(8) :: x, ff
  !
  REAL(8) :: pre1
  INTEGER :: ii
  !
  pre1 = 1.d0/sqrt(N*L)
  ff = 0.d0
  DO ii = 1,N
    ff = ff + cos( PI*( 2*ii - N - 1 )*( x - grid(ibf) )/ L )
  ENDDO
  ff = ff*pre1
END FUNCTION
