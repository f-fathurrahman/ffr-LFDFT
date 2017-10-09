SUBROUTINE dealloc_states()
  USE m_states
  IMPLICIT NONE 
  IF( allocated( KS_evecs ) ) DEALLOCATE(KS_evecs)
  IF( allocated( KS_evals ) ) DEALLOCATE(KS_evals)
  IF( allocated( Focc ) ) DEALLOCATE(Focc)
END SUBROUTINE 
