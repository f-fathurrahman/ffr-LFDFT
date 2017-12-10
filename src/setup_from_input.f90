! FIXME: should be renamed setup_system_from_input ??
SUBROUTINE setup_from_input()
  USE m_input_vars
  IMPLICIT NONE

  CALL setup_atoms()
  CALL setup_LF3d()
  CALL setup_PsPot()

END SUBROUTINE
