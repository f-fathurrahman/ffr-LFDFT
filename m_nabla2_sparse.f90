!! NOTES:
!!
!!   A module to store sparse representation of Laplacian matrix
!!   Naming is based on SparseMatrixCSC type defined in
!!   Julia programming language
!!
!! AUTHOR:
!!   
!!   Fadjar Fathurrahman

MODULE m_nabla2_sparse
  
  IMPLICIT NONE 

  INTEGER :: nabla2_NNZ
  REAL(8), ALLOCATABLE :: nabla2_nzval(:)
  INTEGER, ALLOCATABLE :: nabla2_colptr(:)
  INTEGER, ALLOCATABLE :: nabla2_rowval(:)

END MODULE

