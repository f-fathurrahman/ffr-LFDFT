FUNCTION calc_dir_deriv( N, d, g ) RESULT( dir_deriv )
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: d(N), g(N)
  REAL(8) :: dir_deriv
  !
  REAL(8) :: ddot

  dir_deriv = ddot( N, d, 1, g, 1 )
  dir_deriv = 2.d0*dir_deriv
  
END FUNCTION 
