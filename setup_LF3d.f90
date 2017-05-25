SUBROUTINE setup_LF3d()
  USE m_constants, ONLY : ANG2BOHR
  USE m_input_vars, ONLY : A, B, C, nr1, nr2, nr3
  IMPLICIT NONE 

  CALL init_LF3d_p( (/nr1,nr2,nr3/), (/0.d0,0.d0,0.d0/), (/A,B,C/)*ANG2BOHR )
END SUBROUTINE 

