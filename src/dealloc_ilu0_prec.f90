SUBROUTINE dealloc_ilu0_prec()
  
  USE m_ilu0_prec
  IMPLICIT NONE 

  IF( allocated(alu_ilu0) ) DEALLOCATE( alu_ilu0 )
  IF( allocated(jlu_ilu0) ) DEALLOCATE( jlu_ilu0 )
  IF( allocated(ju_ilu0) ) DEALLOCATE( ju_ilu0 )
  IF( allocated(iw_ilu0) ) DEALLOCATE( iw_ilu0 )

END SUBROUTINE 

