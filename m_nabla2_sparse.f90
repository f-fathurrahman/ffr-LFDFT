
! A module to store sparse representation of Laplacian matrix
MODULE m_nabla2_sparse
  
  IMPLICIT NONE 

  INTEGER :: nabla2_NNZ
  REAL(8), ALLOCATABLE :: nabla2_values(:)
  INTEGER, ALLOCATABLE :: nabla2_column(:)
  INTEGER, ALLOCATABLE :: nabla2_rowIdx(:)

END MODULE

