FUNCTION eval_LF1d_sinc( N, grid, ibf, x ) RESULT(ff)
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: grid(N)
  REAL(8) :: ff, x
  INTEGER :: ibf
  !
  REAL(8) :: h
  
  h = grid(2) - grid(1)
  ff = sin( PI*(x - grid(ibf))/h ) / (PI*(x-grid(ibf))) * h /sqrt(h)
END FUNCTION

