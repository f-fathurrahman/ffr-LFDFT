FUNCTION compute_F( t, x_bar, h ) RESULT( f )
  USE m_constants, ONLY : PI
  IMPLICIT NONE 
  REAL(8) :: t, x_bar, h
  REAL(8) :: f
  !
  REAL(8) :: f_re, f_im

  f = 0.d0

  IF(x_bar < 1.d-30) THEN 
    f = sqrt(h) * erf(PI/(2.d0*h*t))
  ELSE 
    CALL Cwrap_faddeeva( PI/(2.d0*h*t), t*x_bar*h, f_re, f_im )
    f = 0.5d0*exp( -t**2 * x_bar**2 * h**2 ) * h * 2.d0*f_re
  ENDIF 
END FUNCTION 



