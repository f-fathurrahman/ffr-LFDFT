SUBROUTINE dealloc_hamiltonian()
  USE m_hamiltonian
  USE m_xc
  IMPLICIT NONE 
  
  IF( allocated(Rhoe) ) DEALLOCATE( Rhoe )
  IF( allocated(V_ps_loc) ) DEALLOCATE( V_ps_loc )
  IF( allocated(V_Hartree) ) DEALLOCATE( V_Hartree )
  IF( allocated(V_xc) ) DEALLOCATE( V_xc )
  
  IF( allocated(betaNL_psi) ) DEALLOCATE(betaNL_psi)

  IF( allocated(EPS_XC) ) DEALLOCATE( EPS_XC )

END SUBROUTINE 
