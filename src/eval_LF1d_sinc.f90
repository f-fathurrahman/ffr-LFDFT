FUNCTION eval_LF1d_sinc( N, grid, ibf, x ) RESULT(ff)
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: grid(N)
  REAL(8) :: ff, x, dx
  INTEGER :: ibf
  REAL(8), PARAMETER :: SMALL = 1.d-10
  !
  REAL(8) :: h
  
  h = grid(2) - grid(1)
  dx = x - grid(ibf)
  IF( abs(dx) < SMALL ) THEN 
    dx = SMALL
  ENDIF 
  ff = sin( PI*dx/h ) / (PI*dx) * h /sqrt(h)
END FUNCTION

