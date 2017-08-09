SUBROUTINE setup_LF3d()
  USE m_constants, ONLY : ANG2BOHR
  USE m_input_vars, ONLY : A, B, C, nr1, nr2, nr3, assume_isolated
  IMPLICIT NONE 
  REAL(8) :: hh(3)
  
  IF( assume_isolated == 'sinc' ) THEN 
    hh(1) = A/(nr1-1)
    hh(2) = B/(nr2-1)
    hh(3) = C/(nr3-1)
    CALL init_LF3d_sinc( (/nr1,nr2,nr3/), hh*ANG2BOHR )
  ELSE 
    CALL init_LF3d_p( (/nr1,nr2,nr3/), (/0.d0,0.d0,0.d0/), (/A,B,C/)*ANG2BOHR )
  ENDIF 
END SUBROUTINE 

