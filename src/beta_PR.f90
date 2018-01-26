!!> Polak-Ribiere formula
!!>        beta = sum( (g-g_old)*Kg ) / sum( g_old * Kg_old )
FUNCTION beta_PR( N, g, g_old, Kg, Kg_old ) RESULT(beta)
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: g(N), g_old(N), Kg(N), Kg_old(N)
  REAL(8) :: num, denum, beta
  !
  INTEGER :: i
  REAL(8) :: ddot
  
  num = 0.d0
  DO i = 1,N
    num = num + ( g(i) - g_old(i) ) * Kg(i)
  ENDDO 
  denum = ddot( N, g_old, 1, Kg_old, 1 )
  beta = num/denum
END FUNCTION 
