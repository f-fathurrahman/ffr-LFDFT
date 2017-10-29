SUBROUTINE cleanup_ffr_LFDFT()

  USE m_options, ONLY : I_POISSON_SOLVE
  IMPLICIT NONE 

  !
  CALL dealloc_states()

  IF( I_POISSON_SOLVE == 1 ) THEN
    CALL dealloc_Poisson_solve_ISF()
  ELSEIF( I_POISSON_SOLVE == 2 ) THEN
    CALL dealloc_Poisson_solve_DAGE()
  ENDIF

  CALL dealloc_nabla2_sparse()
  CALL dealloc_ilu0_prec()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()
  CALL dealloc_PsPot()
  CALL dealloc_atoms()
END SUBROUTINE 
