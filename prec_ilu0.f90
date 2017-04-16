SUBROUTINE prec_ilu0( v, Kv )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_ilu0_prec

  IMPLICIT NONE 
  REAL(8) :: v(Npoints)
  REAL(8) :: Kv(Npoints)

  CALL lusol( Npoints, v, Kv, alu_ilu0, jlu_ilu0, ju_ilu0 )

END SUBROUTINE 
