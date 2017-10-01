SUBROUTINE dealloc_nabla2_sparse()
  USE m_nabla2_sparse
  IMPLICIT NONE

  IF( allocated(nabla2_nzval) ) DEALLOCATE( nabla2_nzval )
  IF( allocated(nabla2_colptr) ) DEALLOCATE( nabla2_colptr )
  IF( allocated(nabla2_rowval) ) DEALLOCATE( nabla2_rowval )

END SUBROUTINE 
