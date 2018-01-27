FUNCTION beta_DY( N, g, g_old, Kg, d_old ) RESULT( beta )
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: g(N), g_old(N)
  REAL(8) :: Kg(N), d_old(N)
  REAL(8) :: beta
  !
  REAL(8) :: num, denum
  REAL(8) :: ddot
  
  num = ddot( N, g, 1, Kg, 1 )
  denum = ddot( N, g, 1, d_old, 1 ) - ddot( N, g_old, 1, d_old, 1 )
  beta = num/denum
END FUNCTION 

