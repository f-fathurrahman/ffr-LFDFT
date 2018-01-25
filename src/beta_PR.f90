FUNCTION beta_PR( N, g, g_old, Kg, Kg_old ) RESULT(beta)
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: g(N), g_old(N), Kg(N), Kg_old(N)
  REAL(8) :: num, denum, beta
  !
  INTEGER :: i
  
  num = 0.d0
  denum = 0.d0
  DO i = 1,N
    num = num + ( g(i) - g_old(i) ) * Kg(i)
    denum = denum + g_old(i) * Kg_old(i)
  ENDDO 
  beta = num/denum
END FUNCTION 
