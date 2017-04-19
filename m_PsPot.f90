MODULE m_PsPot

  USE m_Ps_HGH, ONLY : Ps_HGH_Params_T
  IMPLICIT NONE 

  CHARACTER(64) :: PsPot_Dir = './HGH/'
  CHARACTER(64), ALLOCATABLE :: PsPot_FilePath(:)

  TYPE(Ps_HGH_Params_T), ALLOCATABLE :: Ps_HGH_Params(:)

END MODULE 
