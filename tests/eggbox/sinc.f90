FUNCTION sinc( x ) RESULT( f )
  USE m_constants, ONLY : PI
  IMPLICIT NONE 
  REAL(8) :: x, f

  f = sin( PI*x )/ (PI*x)
END FUNCTION 

