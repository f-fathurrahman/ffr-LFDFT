FUNCTION interp_LF1d_p( N, L, grid, coefs, x ) RESULT(f)

  IMPLICIT NONE 
  !
  INTEGER :: N
  REAL(8) :: L
  REAL(8) :: grid(N)
  REAL(8) :: coefs(N)
  REAL(8) :: x
  REAL(8) :: f
  !
  INTEGER :: i
  REAL(8) :: eval_LF1d_p

  f = 0.d0
  DO i = 1,N
    f = f + eval_LF1d_p( N, L, grid, i, x )*coefs(i)
  ENDDO 

END FUNCTION 

