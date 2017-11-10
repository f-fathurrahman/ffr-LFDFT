FUNCTION sinc( x ) RESULT( f )
  USE m_constants, ONLY : PI
  IMPLICIT NONE 
  REAL(8) :: SMALL = 1.d-8
  REAL(8) :: x, f

  IF( abs(x) < SMALL ) THEN 
    f = 1.d0
  ELSE 
    f = sin( PI*x )/ (PI*x)
  ENDIF 

END FUNCTION 

