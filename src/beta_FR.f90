!!> Fletcher-Reeves formula
!!>        beta = sum( g * Kg ) / sum( g_old * Kg_old )
FUNCTION beta_FR(N, g, g_old, Kg, Kg_old ) RESULT(beta)
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: g(N)
  REAL(8) :: g_old(N)
  REAL(8) :: Kg(N)
  REAL(8) :: Kg_old(N)
  REAL(8) :: beta
  !
  REAL(8) :: ddot
  REAL(8) :: num, denum

  num = ddot( N, g, 1, Kg, 1 )
  denum = ddot( N, g_old, 1, Kg_old, 1 )
  beta = num/denum

END FUNCTION 
