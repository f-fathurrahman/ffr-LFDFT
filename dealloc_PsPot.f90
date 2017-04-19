SUBROUTINE dealloc_PsPot()

  USE m_PsPot
  IMPLICIT NONE 

  IF( allocated(PsPot_FilePath) ) DEALLOCATE(PsPot_FilePath)
  IF( allocated(Ps_HGH_Params) ) DEALLOCATE(Ps_HGH_Params)

END SUBROUTINE 

