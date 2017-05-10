FUNCTION eval_LF1d_c( N, L, A, grid, ibf, x ) RESULT(ff)
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: L
  REAL(8) :: A
  REAL(8) :: grid(N)
  REAL(8) :: ff
  INTEGER :: ibf
  REAL(8) :: x
  !
  REAL(8) :: ki, pre
  INTEGER :: ii

  pre = 2.d0/sqrt( (N+1)*L )
  ff = 0.d0
  DO ii = 1, N
    ki = PI*ii/L
    ff = ff + sin( ki*(grid(ibf) - A) ) * sin( ki*(x - A) )
  ENDDO
  ff = ff*pre

END FUNCTION

