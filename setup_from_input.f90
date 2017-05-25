SUBROUTINE setup_from_input()
  USE m_input_vars
  IMPLICIT NONE 

  CALL setup_atoms()
  CALL setup_PsPot()
  CALL setup_LF3d()

END SUBROUTINE 

