SUBROUTINE dealloc_nabla2_sparse()
  USE m_nabla2_sparse
  IMPLICIT NONE

  IF( allocated(nabla2_values) ) DEALLOCATE( nabla2_values )
  IF( allocated(nabla2_column) ) DEALLOCATE( nabla2_column )
  IF( allocated(nabla2_rowIdx) ) DEALLOCATE( nabla2_rowIdx )

END SUBROUTINE 
